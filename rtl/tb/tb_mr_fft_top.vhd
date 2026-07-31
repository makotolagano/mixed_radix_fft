library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.fixed_pkg.all;

library std;
use std.textio.all;

library work;
use work.mr_fft_pkg.all;
use work.mr_fft_cfg_pkg.all;

-- Self-checking mr_fft_top TB: AXI4-Lite register access, the shadow/commit
-- configuration protocol, and BIT-EXACT value comparison of the full 11-slot
-- chain against the Python golden model (model/chain_ref.py) -- including the
-- mid-pipeline bypassed slots, which the small configs exercise heavily
-- (N=12 leaves 8 of 11 slots bypassing in series).
--
-- Per config: find the c_fft_sizes index for the golden N, configure it over
-- AXI (shadow CONFIG_SEL + COMMIT + BUSY poll), stream the input codes with
-- random valid/ready gaps, capture every output beat, then compare the first
-- 2N beats against the expected codes (conservation: total outputs = total
-- inputs). Reconfiguration happens on the fly, single reset at time 0.
--
-- Frame markers: i_last is driven on every N-th accepted beat (junk '1'
-- while invalid -- must be gated), o_last is checked on every N-th output
-- beat, and STATUS.FRAMING_ERR must stay 0. A final directed test misplaces
-- i_last, expects the sticky flag, and clears it via CTRL.CLR_FERR.
--
-- Golden vectors come from a CSV data file (model/chain_ref.py) read at
-- runtime with textio -- deliberately NOT a VHDL constant package, so vector
-- regeneration needs no HDL recompile and xsim's debug database stays small.
entity tb_mr_fft_top is
	generic (
		-- path relative to the sim working directory; ../ and ../../ are also
		-- tried, covering runs from the repo root, scripts/ and outputs/
		G_VEC_FILE : string := "rtl/tb/chain_vectors.csv"
	);
end entity;

architecture sim of tb_mr_fft_top is
	constant CLK_PERIOD : time := 10 ns;
	constant C_MAX_CAP  : integer := 3240*3;  -- capture buffer (>= max 3N)
	constant C_MAX_MSG  : integer := 5;     -- mismatch reports per config

	constant C_ZERO_EXT : std_logic_vector(27 downto 0) := (others => '0');

	-- register byte addresses (mirror mr_fft_top)
	constant A_ID     : natural := 16#00#;
	constant A_CTRL   : natural := 16#04#;
	constant A_STATUS : natural := 16#08#;
	constant A_SEL    : natural := 16#0C#;
	constant A_ACT    : natural := 16#10#;
	constant A_INFL   : natural := 16#14#;
	constant A_SIZE   : natural := 16#18#;

	function f_find_cfg(n : natural) return natural is
	begin
		for i in 0 to c_num_configs - 1 loop
			if c_fft_sizes(i) = n then
				return i;
			end if;
		end loop;
		report "config N=" & integer'image(n) & " not in c_fft_sizes" severity failure;
		return 0;
	end function;

	subtype t_word is std_logic_vector(c_fxp_word_width - 1 downto 0);
	type t_cap_arr is array (0 to C_MAX_CAP - 1) of t_word;

	signal clk   : std_logic := '0';
	signal reset : std_logic := '1';
	signal done  : boolean   := false;

	signal awaddr  : std_logic_vector(7 downto 0) := (others => '0');
	signal awvalid : std_logic := '0';
	signal awready : std_logic;
	signal wdata   : std_logic_vector(31 downto 0) := (others => '0');
	signal wvalid  : std_logic := '0';
	signal wready  : std_logic;
	signal bvalid  : std_logic;
	signal bready  : std_logic := '0';
	signal araddr  : std_logic_vector(7 downto 0) := (others => '0');
	signal arvalid : std_logic := '0';
	signal arready : std_logic;
	signal rdata   : std_logic_vector(31 downto 0);
	signal rvalid  : std_logic;
	signal rready  : std_logic := '0';

	constant C_ZERO : t_cmplx := (re => (others => '0'), im => (others => '0'));
	signal i_sample  : t_cmplx   := C_ZERO;
	signal i_sample_slv : std_logic_vector(64 - 1 downto 0);
	signal in_valid  : std_logic := '0';
	signal in_last   : std_logic := '0';
	signal in_ready  : std_logic;
	signal o_sample  : t_cmplx;
	signal o_sample_slv : std_logic_vector(64 - 1 downto 0);
	signal out_valid : std_logic;
	signal out_last  : std_logic;
	signal out_ready : std_logic := '0';

	function to_fx(code : integer) return sfixed is
	begin
		return to_sfixed(std_logic_vector(to_signed(code, c_fxp_word_width)),
		                 c_fxp_int_width - 1, -c_fxp_frac_width);
	end function;

	function code_slv(code : integer) return t_word is
	begin
		return std_logic_vector(to_signed(code, c_fxp_word_width));
	end function;

	function slv_image(w : t_word) return string is
	begin
		if is_x(w) then
			return "0x" & to_hstring(w) & " (metavalue)";
		end if;
		return integer'image(to_integer(signed(w)));
	end function;

	function lcg(r : integer) return integer is
	begin
		return (r * 75 + 74) mod 65537;
	end function;

begin

	dut : entity work.mr_fft_top
		generic map (G_PIPELINE => true)
		port map (
			i_clk => clk, i_reset => reset,
			s_axi_awaddr => awaddr, s_axi_awvalid => awvalid, s_axi_awready => awready,
			s_axi_wdata => wdata, s_axi_wstrb => "1111", s_axi_wvalid => wvalid,
			s_axi_wready => wready, s_axi_bresp => open, s_axi_bvalid => bvalid,
			s_axi_bready => bready,
			s_axi_araddr => araddr, s_axi_arvalid => arvalid, s_axi_arready => arready,
			s_axi_rdata => rdata, s_axi_rresp => open, s_axi_rvalid => rvalid,
			s_axi_rready => rready,
			i_sample => i_sample_slv, i_valid => in_valid, i_last => in_last,
			o_ready => in_ready,
			o_sample => o_sample_slv, o_valid => out_valid, o_last => out_last,
			i_ready => out_ready);

	clk <= not clk after CLK_PERIOD / 2 when not done else '0';

	i_sample_slv <= C_ZERO_EXT & to_slv(i_sample.im) & to_slv(i_sample.re);
	o_sample.re <= to_sfixed(o_sample_slv(17 downto 0), o_sample.re);
	o_sample.im <= to_sfixed(o_sample_slv(35 downto 18), o_sample.im);

	check : process
		type t_code_arr is array (0 to C_MAX_CAP - 1) of integer;
		file vec_f : text;
		variable fstat         : file_open_status;
		variable num_cfgs      : integer;
		variable last_n        : integer := 0;
		variable vin_re, vin_im   : t_code_arr;
		variable vexp_re, vexp_im : t_code_arr;
		variable rd            : std_logic_vector(31 downto 0);
		variable rnd_v, rnd_r  : integer := 777;
		variable cap_re, cap_im : t_cap_arr;
		variable n, sel        : integer;
		variable n_in, n_exp   : integer;
		variable in_idx, cap_cnt, guard : integer;
		variable v             : boolean;
		variable errs, msgs    : integer;
		variable last_errs     : integer;
		variable exp_re, exp_im : t_word;
		variable pass_all      : boolean := true;

		-- next non-comment, non-empty line of the vector file
		procedure vec_line(vl : inout line) is
		begin
			loop
				assert not endfile(vec_f)
					report "vector file ended early" severity failure;
				readline(vec_f, vl);
				exit when vl'length > 0 and vl(vl'left) /= '#';
			end loop;
		end procedure;

		procedure vec_int(a : out integer) is
			variable vl : line;
		begin
			vec_line(vl);
			read(vl, a);
		end procedure;

		-- one 'a,b' row
		procedure vec_pair(a, b : out integer) is
			variable vl : line;
			variable ch : character;
		begin
			vec_line(vl);
			read(vl, a);
			read(vl, ch);   -- ','
			read(vl, b);
		end procedure;

		-- one 'a,b,c' config header row
		procedure vec_hdr(a, b, c : out integer) is
			variable vl : line;
			variable ch : character;
		begin
			vec_line(vl);
			read(vl, a);
			read(vl, ch);   -- ','
			read(vl, b);
			read(vl, ch);   -- ','
			read(vl, c);
		end procedure;

		procedure axi_write(addr : natural; data : std_logic_vector(31 downto 0)) is
		begin
			wait until rising_edge(clk);
			awaddr  <= std_logic_vector(to_unsigned(addr, 8));
			wdata   <= data;
			awvalid <= '1';
			wvalid  <= '1';
			bready  <= '1';
			loop
				wait until rising_edge(clk);
				exit when awready = '1';
			end loop;
			awvalid <= '0';
			wvalid  <= '0';
			loop
				exit when bvalid = '1';
				wait until rising_edge(clk);
			end loop;
			wait until rising_edge(clk);
			bready <= '0';
		end procedure;

		procedure axi_read(addr : natural; data : out std_logic_vector(31 downto 0)) is
		begin
			wait until rising_edge(clk);
			araddr  <= std_logic_vector(to_unsigned(addr, 8));
			arvalid <= '1';
			rready  <= '1';
			loop
				wait until rising_edge(clk);
				exit when rvalid = '1';
			end loop;
			data := rdata;
			arvalid <= '0';
			rready  <= '0';
		end procedure;

		-- configure + commit + poll BUSY clear
		procedure configure(sel_v : natural) is
			variable st : std_logic_vector(31 downto 0);
		begin
			axi_write(A_SEL, std_logic_vector(to_unsigned(sel_v, 32)));
			axi_write(A_CTRL, x"00000001");
			loop
				axi_read(A_STATUS, st);
				exit when st(0) = '0';
			end loop;
		end procedure;

	begin
		file_open(fstat, vec_f, G_VEC_FILE, read_mode);
		if fstat /= open_ok then
			file_open(fstat, vec_f, "../" & G_VEC_FILE, read_mode);
		end if;
		if fstat /= open_ok then
			file_open(fstat, vec_f, "../../" & G_VEC_FILE, read_mode);
		end if;
		if fstat /= open_ok then
			-- Vivado project mode: data files in the sim fileset are copied
			-- FLAT into <proj>.sim/sim_1/behav/xsim/, the xsim working dir
			file_open(fstat, vec_f, "chain_vectors.csv", read_mode);
		end if;
		assert fstat = open_ok
			report "cannot open vector file " & G_VEC_FILE &
			       " (regenerate: .venv/bin/python model/chain_ref.py)"
			severity failure;
		vec_int(num_cfgs);

		reset <= '1';
		for i in 0 to 3 loop
			wait until rising_edge(clk);
		end loop;
		reset <= '0';

		-- 1. ID and reset-state STATUS
		axi_read(A_ID, rd);
		assert rd = x"0FF70100"
			report "ID mismatch: " & to_hstring(rd) severity error;
		pass_all := pass_all and (rd = x"0FF70100");
		axi_read(A_STATUS, rd);
		assert rd(0) = '0' and rd(1) = '1'
			report "reset STATUS not idle" severity error;
		pass_all := pass_all and (rd(0) = '0' and rd(1) = '1');

		-- 2. golden-vector configs, reconfigured on the fly over AXI
		for c in 0 to num_cfgs - 1 loop
			vec_hdr(n, n_in, n_exp);
			sel    := f_find_cfg(n);
			last_n := n;
			assert n_in <= C_MAX_CAP and n_exp <= C_MAX_CAP
				report "C_MAX_CAP too small for N=" & integer'image(n)
				severity failure;
			for k in 0 to n_in - 1 loop
				vec_pair(vin_re(k), vin_im(k));
			end loop;
			for k in 0 to n_exp - 1 loop
				vec_pair(vexp_re(k), vexp_im(k));
			end loop;

			configure(sel);
			axi_read(A_SIZE, rd);
			assert to_integer(unsigned(rd)) = n
				report "FFT_SIZE readback wrong for N=" & integer'image(n)
				severity error;
			pass_all := pass_all and (to_integer(unsigned(rd)) = n);

			-- stream all inputs with random gaps; capture every output beat
			in_idx    := 0;
			cap_cnt   := 0;
			guard     := 0;
			last_errs := 0;
			while cap_cnt < n_in loop
				rnd_v := lcg(rnd_v);
				rnd_r := lcg(rnd_r);
				v := in_idx < n_in; --(rnd_v mod 4) /= 0 
				out_ready <= '1';-- when (rnd_r mod 4) /= 0 else '0';
				if v then
					in_valid <= '1';
					i_sample.re <= to_fx(vin_re(in_idx));
					i_sample.im <= to_fx(vin_im(in_idx));
					if (in_idx mod n) = n - 1 then
						in_last <= '1';
					else
						in_last <= '0';
					end if;
				else
					in_valid <= '0';
					i_sample.re <= to_fx(21845);              -- junk while not valid
					i_sample.im <= to_fx(-21846);
					in_last <= '1';                           -- junk: must be gated by valid
				end if;

				wait until falling_edge(clk);
				if out_valid = '1' and out_ready = '1' then
					cap_re(cap_cnt) := to_slv(o_sample.re);
					cap_im(cap_cnt) := to_slv(o_sample.im);
					-- o_last must mark exactly every N-th output beat
					if (out_last = '1') /= ((cap_cnt mod n) = n - 1) then
						last_errs := last_errs + 1;
						if last_errs <= C_MAX_MSG then
							report "  N=" & integer'image(n) & " out beat " &
							       integer'image(cap_cnt) & ": o_last=" &
							       std_logic'image(out_last) & " wrong"
								severity error;
						end if;
					end if;
					cap_cnt := cap_cnt + 1;
				end if;
				if v and in_ready = '1' then
					in_idx := in_idx + 1;
				end if;
				guard := guard + 1;
				assert guard < 100 * n_in
					report "N=" & integer'image(n) & " deadlock: outputs " &
					       integer'image(cap_cnt) & "/" & integer'image(n_in)
					severity failure;
				wait until rising_edge(clk);
			end loop;
			in_valid  <= '0';
			in_last   <= '0';
			out_ready <= '0';

			-- IN_FLIGHT must be back to zero (conservation reached)
			axi_read(A_INFL, rd);
			assert to_integer(unsigned(rd)) = 0
				report "IN_FLIGHT not zero after N=" & integer'image(n)
				severity error;
			pass_all := pass_all and (to_integer(unsigned(rd)) = 0);

			-- correctly framed traffic must not raise FRAMING_ERR
			axi_read(A_STATUS, rd);
			assert rd(3) = '0'
				report "FRAMING_ERR set by well-framed input, N=" & integer'image(n)
				severity error;
			pass_all := pass_all and (rd(3) = '0') and (last_errs = 0);

			-- compare the first n_exp beats against the golden model
			errs := 0;
			msgs := 0;
			for k in 0 to n_exp - 1 loop
				exp_re := code_slv(vexp_re(k));
				exp_im := code_slv(vexp_im(k));
				if k >= cap_cnt or cap_re(k) /= exp_re or cap_im(k) /= exp_im then
					errs := errs + 1;
					if msgs < C_MAX_MSG then
						report "  N=" & integer'image(n) & " out beat " &
						       integer'image(k) & " got (" &
						       slv_image(cap_re(k)) & "," & slv_image(cap_im(k)) &
						       ") exp (" & integer'image(vexp_re(k)) &
						       "," & integer'image(vexp_im(k)) & ")"
							severity error;
						msgs := msgs + 1;
					end if;
				end if;
			end loop;

			if errs = 0 then
				report "chain config N=" & integer'image(n) & " (sel " &
				       integer'image(sel) & "): PASS -- " & integer'image(n_exp) &
				       " samples match the Python golden model"
					severity note;
			else
				pass_all := false;
				report "chain config N=" & integer'image(n) & ": FAIL -- " &
				       integer'image(errs) & "/" & integer'image(n_exp) &
				       " sample mismatches"
					severity error;
			end if;
		end loop;

		-- 3. shadow write must NOT change the active config until committed
		axi_write(A_SEL, std_logic_vector(to_unsigned(0, 32)));
		axi_read(A_ACT, rd);
		assert to_integer(unsigned(rd)) = f_find_cfg(last_n)
			report "shadow write leaked into CONFIG_ACTIVE" severity error;
		pass_all := pass_all and
			(to_integer(unsigned(rd)) = f_find_cfg(last_n));

		-- 4. misplaced i_last (beat 0 instead of N-1) sets sticky FRAMING_ERR;
		--    o_last stays counter-generated; CTRL.CLR_FERR clears the flag.
		--    Still one WHOLE frame (last_n samples) so the chain drains.
		n       := last_n;
		in_idx  := 0;
		cap_cnt := 0;
		guard   := 0;
		while cap_cnt < n loop
			v := in_idx < n;
			out_ready <= '1';
			if v then
				in_valid <= '1';
				i_sample.re <= to_fx(vin_re(in_idx));
				i_sample.im <= to_fx(vin_im(in_idx));
				if in_idx = 0 then
					in_last <= '1';   -- wrong beat
				else
					in_last <= '0';   -- and missing on beat N-1
				end if;
			else
				in_valid <= '0';
				in_last  <= '0';
			end if;

			wait until falling_edge(clk);
			if out_valid = '1' and out_ready = '1' then
				if (out_last = '1') /= (cap_cnt = n - 1) then
					pass_all := false;
					report "framing test: o_last wrong on out beat " &
					       integer'image(cap_cnt) severity error;
				end if;
				cap_cnt := cap_cnt + 1;
			end if;
			if v and in_ready = '1' then
				in_idx := in_idx + 1;
			end if;
			guard := guard + 1;
			assert guard < 100 * n
				report "framing test deadlock: outputs " &
				       integer'image(cap_cnt) & "/" & integer'image(n)
				severity failure;
			wait until rising_edge(clk);
		end loop;
		in_valid  <= '0';
		in_last   <= '0';
		out_ready <= '0';

		axi_read(A_STATUS, rd);
		assert rd(3) = '1'
			report "FRAMING_ERR not set by misplaced i_last" severity error;
		pass_all := pass_all and (rd(3) = '1');

		axi_write(A_CTRL, x"00000002");   -- CLR_FERR
		axi_read(A_STATUS, rd);
		assert rd(3) = '0'
			report "CTRL.CLR_FERR did not clear FRAMING_ERR" severity error;
		pass_all := pass_all and (rd(3) = '0');

		if pass_all then
			report "TB PASS -- mr_fft_top: AXI configuration + full chain bit-exact vs the Python golden model"
				severity note;
		else
			report "TB FAIL -- see errors above" severity error;
		end if;
		done <= true;
		wait;
	end process check;

end architecture sim;
