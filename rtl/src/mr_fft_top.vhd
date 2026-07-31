library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.fixed_pkg.all;

library work;
use work.mr_fft_pkg.all;
use work.mr_fft_cfg_pkg.all;

-- ---------------------------------------------------------------------------
-- Mixed-radix FFT top: the mid-bypass chain plus an AXI4-Lite configuration
-- interface. Single clock domain (AXI and stream share i_clk).
--
-- Register map (32-bit registers, byte addresses):
--   0x00  ID           RO  x"0FF70100" (FFT core, v1.0)
--   0x04  CTRL         WO  bit0 COMMIT (self-clearing), bit1 CLR_FERR
--   0x08  STATUS       RO  bit0 BUSY (commit pending), bit1 IDLE (in_flight=0),
--                          bit3 FRAMING_ERR (sticky, clear via CTRL.CLR_FERR)
--   0x0C  CONFIG_SEL   RW  SHADOW config select (index into c_fft_sizes)
--   0x10  CONFIG_ACTIVE RO active (committed) config select
--   0x14  IN_FLIGHT    RO  input beats minus output beats inside the chain
--   0x18  FFT_SIZE     RO  N of the ACTIVE config (from c_fft_sizes)
--
-- Configuration protocol (encodes the drained-switch contract in hardware):
-- writes to CONFIG_SEL land in a shadow register only. Writing CTRL.COMMIT
-- gates the input stream (o_ready forced low), waits until the chain is
-- drained (IN_FLIGHT = 0 -- exact, by the chain's beat conservation), then
-- transfers shadow -> active in one cycle and reopens the input. Software:
-- write CONFIG_SEL, write COMMIT, poll STATUS.BUSY = 0 (or just keep
-- streaming: the gate handles the boundary). NOTE: commit only after WHOLE
-- frames have been offered -- a partial frame parks samples in the delay
-- FIFOs and IN_FLIGHT never reaches zero (STATUS makes this visible).
--
-- Frame marker (AXI-Stream tlast semantics): o_last is REGENERATED from a
-- mod-N output beat counter (N of the ACTIVE config), never forwarded from
-- the input -- exact because a commit only lands with the chain drained on a
-- frame boundary. i_last cannot steer the datapath (the delay FIFOs are
-- configured for N), so it is only checked: an accepted beat where i_last
-- disagrees with the input beat counter sets the sticky STATUS.FRAMING_ERR.
-- ---------------------------------------------------------------------------
entity mr_fft_top is
	generic (
		G_PIPELINE : boolean := true
	);
	port (
		i_clk   : in  std_logic;
		i_reset : in  std_logic;

		-- AXI4-Lite slave (configuration)
		s_axi_awaddr  : in  std_logic_vector(7 downto 0);
		s_axi_awvalid : in  std_logic;
		s_axi_awready : out std_logic;
		s_axi_wdata   : in  std_logic_vector(31 downto 0);
		s_axi_wstrb   : in  std_logic_vector(3 downto 0);
		s_axi_wvalid  : in  std_logic;
		s_axi_wready  : out std_logic;
		s_axi_bresp   : out std_logic_vector(1 downto 0);
		s_axi_bvalid  : out std_logic;
		s_axi_bready  : in  std_logic;
		s_axi_araddr  : in  std_logic_vector(7 downto 0);
		s_axi_arvalid : in  std_logic;
		s_axi_arready : out std_logic;
		s_axi_rdata   : out std_logic_vector(31 downto 0);
		s_axi_rresp   : out std_logic_vector(1 downto 0);
		s_axi_rvalid  : out std_logic;
		s_axi_rready  : in  std_logic;

		-- input stream handshake
		i_sample : in  std_logic_vector(63 downto 0);
		i_valid  : in  std_logic;
		i_last   : in  std_logic;
		o_ready  : out std_logic;

		-- output stream handshake
		o_sample : out std_logic_vector(63 downto 0);
		o_valid  : out std_logic;
		o_last   : out std_logic;
		i_ready  : in  std_logic
	);
end entity mr_fft_top;

architecture rtl of mr_fft_top is

	constant C_ID : std_logic_vector(31 downto 0) := x"0FF70100";

	constant C_SEL_W : natural := clogb2(c_num_configs);

	-- register byte addresses (word aligned)
	constant C_ADDR_ID         : natural := 16#00#;
	constant C_ADDR_CTRL       : natural := 16#04#;
	constant C_ADDR_STATUS     : natural := 16#08#;
	constant C_ADDR_CONFIG_SEL : natural := 16#0C#;
	constant C_ADDR_CONFIG_ACT : natural := 16#10#;
	constant C_ADDR_IN_FLIGHT  : natural := 16#14#;
	constant C_ADDR_FFT_SIZE   : natural := 16#18#;
	
	constant C_ZERO_EXT : std_logic_vector(27 downto 0) := (others => '0');

	-- generous upper bound on words inside the chain (FIFOs + skids + pipes)
	constant C_INFLIGHT_W : natural := 16;

	-- frame beat counters count 0 .. N-1 (N <= c_max_fft_size)
	constant C_BEAT_W : natural := clogb2(c_max_fft_size);

	-- configuration registers
	signal shadow_sel : unsigned(C_SEL_W - 1 downto 0);
	signal active_sel : unsigned(C_SEL_W - 1 downto 0);
	signal commit_pnd : std_logic;   -- commit pending: input gated, draining
	signal fft_ifft_shadow_sel : std_logic; -- FFT/IFFT select
	signal fft_ifft_active_sel : std_logic; -- FFT/IFFT select

	-- stream bookkeeping
	signal in_flight   : unsigned(C_INFLIGHT_W - 1 downto 0);
	signal chain_in_valid : std_logic;
	signal chain_ready : std_logic;
	signal chain_valid : std_logic;
	signal top_ready   : std_logic;
	signal in_beat     : std_logic;
	signal out_beat    : std_logic;

	-- frame position tracking (tlast generation / input framing check)
	signal in_cnt      : unsigned(C_BEAT_W - 1 downto 0);
	signal out_cnt     : unsigned(C_BEAT_W - 1 downto 0);
	signal active_n_m1 : unsigned(C_BEAT_W - 1 downto 0);   -- N-1, active config
	signal framing_err : std_logic;

	-- AXI4-Lite slave (single outstanding transaction)
	signal awready, wready, bvalid   : std_logic;
	signal arready, rvalid           : std_logic;
	signal rdata                     : std_logic_vector(31 downto 0);
	signal wr_beat                   : std_logic;
	
	signal in_sample : t_cmplx;
	signal out_sample : t_cmplx;

begin

	in_sample.re <= to_sfixed(i_sample(17 downto 0), in_sample.re);
	in_sample.im <= to_sfixed(i_sample(35 downto 18), in_sample.im);
	
	o_sample <= C_ZERO_EXT & to_slv(out_sample.im) & to_slv(out_sample.re);

	-- ------------------------------------------------------------------
	-- chain + stream gating
	-- ------------------------------------------------------------------
	CHAIN_INST: entity work.mr_fft_chain
		generic map (G_PIPELINE => G_PIPELINE)
		port map (
			i_clk    => i_clk,
			i_reset  => i_reset,
			i_config_sel => std_logic_vector(active_sel),
			i_sample => in_sample,
			i_valid  => chain_in_valid,
			o_ready  => chain_ready,
			o_sample => out_sample,
			o_valid  => chain_valid,
			i_ready  => i_ready
		);

	-- input gate: a pending commit refuses new samples so the chain drains
	chain_in_valid <= i_valid and not commit_pnd;
	top_ready <= chain_ready and not commit_pnd;
	o_ready   <= top_ready;
	o_valid   <= chain_valid;

	in_beat  <= i_valid and top_ready;
	out_beat <= chain_valid and i_ready;

	PROC_IN_FLIGHT: process(i_clk)
	begin
		if rising_edge(i_clk) then
			if i_reset = '1' then
				in_flight <= (others => '0');
			elsif in_beat = '1' and out_beat = '0' then
				in_flight <= in_flight + 1;
			elsif in_beat = '0' and out_beat = '1' then
				in_flight <= in_flight - 1;
			end if;
		end if;
	end process PROC_IN_FLIGHT;

	-- ------------------------------------------------------------------
	-- frame markers: o_last regenerated from a mod-N output beat counter;
	-- i_last only CHECKED against the input beat counter (sticky error)
	-- ------------------------------------------------------------------
	active_n_m1 <= to_unsigned(
	    c_fft_sizes(minimum(to_integer(active_sel), c_num_configs - 1)) - 1,
	    C_BEAT_W);

	o_last <= '1' when out_cnt = active_n_m1 else '0';

	PROC_FRAME: process(i_clk)
	begin
		if rising_edge(i_clk) then
			if i_reset = '1' then
				in_cnt      <= (others => '0');
				out_cnt     <= (others => '0');
				framing_err <= '0';
			else
				-- CTRL.CLR_FERR; a same-cycle framing error below wins
				if wr_beat = '1' and s_axi_wdata(1) = '1'
				   and to_integer(unsigned(s_axi_awaddr)) = C_ADDR_CTRL then
					framing_err <= '0';
				end if;

				if in_beat = '1' then
					if in_cnt = active_n_m1 then
						in_cnt <= (others => '0');
					else
						in_cnt <= in_cnt + 1;
					end if;
					-- i_last must mark exactly the final beat of each frame
					if (i_last = '1') /= (in_cnt = active_n_m1) then
						framing_err <= '1';
					end if;
				end if;

				if out_beat = '1' then
					if out_cnt = active_n_m1 then
						out_cnt <= (others => '0');
					else
						out_cnt <= out_cnt + 1;
					end if;
				end if;

				-- defensive: both counters are provably 0 when a commit
				-- lands (whole-frame drain); clear anyway so a framing slip
				-- cannot survive a reconfiguration
				if commit_pnd = '1' and in_flight = 0 and in_beat = '0' then
					in_cnt  <= (others => '0');
					out_cnt <= (others => '0');
				end if;
			end if;
		end if;
	end process PROC_FRAME;

	-- ------------------------------------------------------------------
	-- commit FSM: gate input -> wait drained -> apply shadow
	-- ------------------------------------------------------------------
	PROC_COMMIT: process(i_clk)
	begin
		if rising_edge(i_clk) then
			if i_reset = '1' then
				shadow_sel <= (others => '0');
				active_sel <= (others => '0');
				commit_pnd <= '0';
				fft_ifft_shadow_sel <= '0';
				fft_ifft_active_sel <= '0';
			else
				-- shadow write (CONFIG_SEL)
				if wr_beat = '1' and to_integer(unsigned(s_axi_awaddr)) = C_ADDR_CONFIG_SEL then
					shadow_sel <= resize(unsigned(s_axi_wdata(C_SEL_W - 1 downto 0)), C_SEL_W);
				end if;

				-- commit request (CTRL.COMMIT)
				if wr_beat = '1' and to_integer(unsigned(s_axi_awaddr)) = C_ADDR_CTRL then
					commit_pnd <= s_axi_wdata(0);
					fft_ifft_shadow_sel <= s_axi_wdata(2);
				end if;

				if commit_pnd = '1' and in_flight = 0 and in_beat = '0' then
					active_sel <= shadow_sel;
					fft_ifft_active_sel <= fft_ifft_shadow_sel;
					commit_pnd <= '0';
				end if;
			end if;
		end if;
	end process PROC_COMMIT;

	-- ------------------------------------------------------------------
	-- AXI4-Lite slave: aw+w accepted together, single outstanding
	-- ------------------------------------------------------------------
	wr_beat <= awready and s_axi_awvalid and s_axi_wvalid;

	PROC_AXI_WR: process(i_clk)
	begin
		if rising_edge(i_clk) then
			if i_reset = '1' then
				awready <= '0';
				wready  <= '0';
				bvalid  <= '0';
			else
				awready <= '0';
				wready  <= '0';
				if bvalid = '1' then
					if s_axi_bready = '1' then
						bvalid <= '0';
					end if;
				elsif awready = '0' and s_axi_awvalid = '1' and s_axi_wvalid = '1' then
					awready <= '1';
					wready  <= '1';
					bvalid  <= '1';
				end if;
			end if;
		end if;
	end process PROC_AXI_WR;

	-- accept aw and w in the same beat (awready/wready pulse together);
	-- register writes happen on wr_beat inside PROC_COMMIT
	s_axi_awready <= awready;
	s_axi_wready  <= wready;
	s_axi_bvalid  <= bvalid;
	s_axi_bresp   <= "00";   -- always OKAY (unknown addresses write nothing)

	PROC_AXI_RD: process(i_clk)
		variable v_addr : natural;
	begin
		if rising_edge(i_clk) then
			if i_reset = '1' then
				arready <= '0';
				rvalid  <= '0';
				rdata   <= (others => '0');
			else
				arready <= '0';
				if rvalid = '1' then
					if s_axi_rready = '1' then
						rvalid <= '0';
					end if;
				elsif arready = '0' and s_axi_arvalid = '1' then
					arready <= '1';
					rvalid  <= '1';
					v_addr  := to_integer(unsigned(s_axi_araddr));
					rdata   <= (others => '0');
					case v_addr is
						when C_ADDR_ID =>
							rdata <= C_ID;
						when C_ADDR_STATUS =>
							rdata(0) <= commit_pnd;
							if in_flight = 0 then
								rdata(1) <= '1';
							end if;
							rdata(2) <= fft_ifft_active_sel;
							rdata(3) <= framing_err;
						when C_ADDR_CONFIG_SEL =>
							rdata(C_SEL_W - 1 downto 0) <= std_logic_vector(shadow_sel);
						when C_ADDR_CONFIG_ACT =>
							rdata(C_SEL_W - 1 downto 0) <= std_logic_vector(active_sel);
						when C_ADDR_IN_FLIGHT =>
							rdata(C_INFLIGHT_W - 1 downto 0) <= std_logic_vector(in_flight);
						when C_ADDR_FFT_SIZE =>
							rdata <= std_logic_vector(to_unsigned(
							    c_fft_sizes(minimum(to_integer(active_sel), c_num_configs - 1)), 32));
						when others =>
							rdata <= (others => '0');
					end case;
				end if;
			end if;
		end if;
	end process PROC_AXI_RD;

	s_axi_arready <= arready;
	s_axi_rvalid  <= rvalid;
	s_axi_rdata   <= rdata;
	s_axi_rresp   <= "00";

end architecture rtl;
