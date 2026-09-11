library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.fixed_pkg.all;
use ieee.fixed_float_types.all;

library work;
use work.mr_fft_pkg.all;

-- Final scaler: o_sample = i_sample * c, c = 2**total_shift / N of the active
-- config, so the output follows the X/N convention. c is real, so one 25x18
-- multiply per component: constant on the DSP A port, data on the B port,
-- rounding constant on the C port.
--
-- pipeline: A operand register, B biased product (M/P regs), C truncate.
-- a credit gated skid keeps o_ready registered.
-- i_scale is quasi-static, the top changes it together with the config.
entity mr_fft_scaler is
	port (
		i_clk   : in  std_logic;
		i_reset : in  std_logic;

		-- scale code of the active config
		i_scale : in  std_logic_vector(c_scale_int + c_scale_frac - 1 downto 0);

		-- input stream handshake
		i_sample : in  t_cmplx;
		i_valid  : in  std_logic;
		o_ready  : out std_logic;

		-- output stream handshake
		o_sample : out t_cmplx;
		o_valid  : out std_logic;
		i_ready  : in  std_logic
	);
end entity mr_fft_scaler;

architecture rtl of mr_fft_scaler is

	constant C_LAT        : natural := 3;
	constant C_SKID_DEPTH : natural := C_LAT + 3;

	subtype t_scale is sfixed(c_scale_int - 1 downto -c_scale_frac);
	-- data x scale product
	subtype t_prod is sfixed(c_fxp_int_width + c_scale_int - 1
	                         downto -(c_fxp_frac_width + c_scale_frac));

	-- half-up rounding constant, half LSB of the data word
	constant c_half : sfixed(0 downto -(c_fxp_frac_width + 1)) :=
		to_sfixed(2.0 ** (-(c_fxp_frac_width + 1)), 0, -(c_fxp_frac_width + 1));

	signal scale_fx : t_scale;

	-- stage A: operand register
	signal a_r : t_cmplx;
	-- stage B: biased products (DSP M/P registers)
	signal p_re_r, p_im_r : t_prod;
	attribute use_dsp : string;
	attribute use_dsp of p_re_r, p_im_r : signal is "yes";

	-- valid sideband (A, B, C)
	signal v : std_logic_vector(1 to C_LAT);

	-- output skid and credit
	signal skid_we, skid_re  : std_logic;
	signal skid_data         : t_cmplx;
	signal skid_out          : t_cmplx;
	signal skid_valid        : std_logic;
	signal out_credit        : integer range 0 to C_SKID_DEPTH;
	signal credit_ok         : std_logic;
	signal launch            : std_logic;

begin

	scale_fx <= to_sfixed(i_scale, c_scale_int - 1, -c_scale_frac);

	credit_ok <= '1' when out_credit < C_SKID_DEPTH else '0';
	o_ready   <= credit_ok;
	launch    <= i_valid and credit_ok;

	PROC_PIPE: process(i_clk)
	begin
		if rising_edge(i_clk) then
			-- stage A
			a_r <= i_sample;

			-- stage B: one real multiply per component, rounding bias on the post-adder
			p_re_r <= resize(a_r.re * scale_fx + c_half, c_fxp_int_width + c_scale_int - 1,
			                 -(c_fxp_frac_width + c_scale_frac), fixed_wrap, fixed_truncate);
			p_im_r <= resize(a_r.im * scale_fx + c_half, c_fxp_int_width + c_scale_int - 1,
			                 -(c_fxp_frac_width + c_scale_frac), fixed_wrap, fixed_truncate);

			-- stage C: truncate to the data word
			skid_data.re <= resize(p_re_r, skid_data.re, fixed_wrap, fixed_truncate);
			skid_data.im <= resize(p_im_r, skid_data.im, fixed_wrap, fixed_truncate);

			-- valid alongside
			if i_reset = '1' then
				v <= (others => '0');
			else
				v <= launch & v(1 to C_LAT - 1);
			end if;
		end if;
	end process PROC_PIPE;

	skid_we <= v(C_LAT);

	SKID_INST: entity work.mr_fft_fifo
		generic map (G_DEPTH => C_SKID_DEPTH)
		port map (
			i_clk       => i_clk,
			i_reset     => i_reset,
			i_wr_en     => skid_we,
			i_wr_sample => skid_data,
			i_rd_en     => skid_re,
			o_rd_sample => skid_out,
			o_rd_valid  => skid_valid,
			o_full      => open
		);

	o_valid  <= skid_valid;
	o_sample <= skid_out;
	skid_re  <= skid_valid and i_ready;

	PROC_CREDIT: process(i_clk)
		variable v_inc, v_dec : boolean;
	begin
		if rising_edge(i_clk) then
			if i_reset = '1' then
				out_credit <= 0;
			else
				v_inc := launch = '1';
				v_dec := skid_valid = '1' and i_ready = '1';
				if v_inc and not v_dec then
					out_credit <= out_credit + 1;
				elsif v_dec and not v_inc then
					out_credit <= out_credit - 1;
				end if;
			end if;
		end if;
	end process PROC_CREDIT;

end architecture rtl;
