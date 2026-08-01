library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.fixed_pkg.all;
use ieee.fixed_float_types.all;

library work;
use work.mr_fft_pkg.all;

-- Final scaler: o_sample = i_sample * c, where c = 2**total_shift / N is the
-- ACTIVE config's residual gain (u2.<c_scale_frac> code on i_scale, from the
-- c_fft_scales elaboration table, s<c_scale_int>.<c_scale_frac>) --
-- turning the chain's X * 2**(-total_shift)
-- output into the classical X/N convention. c is real, so each component is
-- ONE 25x18 multiply: the 25-bit constant rides the DSP A port, the data
-- word the 18-bit B port, and the half-up exit rounding constant (+half LSB
-- of the memory word) rides the post-adder C port, so the M/P registers
-- absorb the pipeline registers.
--
-- Pipeline (free-running, valid sideband, same recipe as mr_fft_rotator):
--   stage A: operand register
--   stage B: biased product (data * c + half) -> DSP M/P register
--   stage C: truncate to the memory word (the ONE rounding, wrap)
-- A credit-gated output skid keeps the handshake registered: o_ready is the
-- (registered) credit state, never a combinational path from i_ready.
--
-- i_scale is quasi-static: change it only when the scaler is drained (the
-- top switches it together with the config at a committed drain boundary).
entity mr_fft_scaler is
	port (
		i_clk   : in  std_logic;
		i_reset : in  std_logic;

		-- s<c_scale_int>.<c_scale_frac> scale code of the active config (quasi-static)
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
	-- data (s2.16) x scale (s4.21) product, plus the rounding bias headroom
	subtype t_prod is sfixed(c_fxp_int_width + c_scale_int - 1
	                         downto -(c_fxp_frac_width + c_scale_frac));

	-- half-up exit rounding constant: +half LSB of the memory word
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

	-- output skid + credit (in-flight words counted -> no overflow possible)
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

			-- stage B: one real multiply per component, trim bias on the
			-- DSP post-adder
			p_re_r <= resize(a_r.re * scale_fx + c_half, c_fxp_int_width + c_scale_int - 1,
			                 -(c_fxp_frac_width + c_scale_frac), fixed_wrap, fixed_truncate);
			p_im_r <= resize(a_r.im * scale_fx + c_half, c_fxp_int_width + c_scale_int - 1,
			                 -(c_fxp_frac_width + c_scale_frac), fixed_wrap, fixed_truncate);

			-- stage C: the single truncation to the memory word (wrap)
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
