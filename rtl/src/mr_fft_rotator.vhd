library ieee;
use ieee.std_logic_1164.all;
use ieee.fixed_pkg.all;
use ieee.fixed_float_types.all;

-- Custom packages
library work;
use work.mr_fft_pkg.all;

-- Rotator: o_sample = i_sample * i_twiddle, bit-exact to the mr_fft_pkg "*"
-- operator (full-precision 4-multiplier complex product, ONE final half-up
-- exit rounding: +half LSB then truncate, wrap).
--
-- G_PIPELINE = false: combinational (original behavior), latency 0.
-- G_PIPELINE = true : rotator_latency(true) = 5 register stages, placed so
-- Vivado absorbs them into the four DSP48E1s:
--   stage A1/A2: operand registers, both   -> DSP A/B input registers
--   stage B: the four 18x18 products       -> DSP M registers
--   stage C: full-precision sum/difference -> DSP post-adders / P registers
--   stage D: rounded resize to the data word (fabric)
-- The datapath free-runs (no enables); the valid bit travels alongside, so
-- the consumer qualifies outputs exactly like the preadder pipeline.
entity mr_fft_rotator is
	generic (
		G_PIPELINE : boolean := true
	);
	port (
		i_clk 		: in  std_logic := '0';
		i_reset 	: in  std_logic := '0';

		i_valid 	: in  std_logic := '1';
		i_sample 	: in  t_cmplx;
		i_twiddle : in  t_cmplx_twiddle;

		o_valid 	: out std_logic;
		o_sample 	: out t_cmplx
	);
end entity mr_fft_rotator;

architecture rtl of mr_fft_rotator is

	-- one 18x18 product: (data int + twiddle int) ints, summed frac bits
	subtype t_prod is sfixed(c_fxp_int_width + c_twiddle_int_width - 1
	                         downto -(c_fxp_frac_width + c_twiddle_frac_width));
	-- product sum/difference: one growth bit (same shape as t_cmplx_mult)
	subtype t_prod_sum is sfixed(c_fxp_int_width + c_twiddle_int_width
	                             downto -(c_fxp_frac_width + c_twiddle_frac_width));
	-- half-up exit rounding constant: +half LSB of the memory word
	constant c_half : sfixed(0 downto -(c_fxp_frac_width + 1)) :=
		to_sfixed(2.0 ** (-(c_fxp_frac_width + 1)), 0, -(c_fxp_frac_width + 1));

begin

	GEN_COMB: if not G_PIPELINE generate
		o_sample <= i_sample * i_twiddle;
		o_valid  <= i_valid;
	end generate GEN_COMB;

	GEN_PIPE: if G_PIPELINE generate
		-- stages A1/A2: operand registers
		signal a_r1, b_r1, a_r : sfixed(c_fxp_int_width - 1 downto -c_fxp_frac_width);
		signal b_r : sfixed(c_fxp_int_width - 1 downto -c_fxp_frac_width);
		signal c_r1, d_r1, c_r : sfixed(c_twiddle_int_width - 1 downto -c_twiddle_frac_width);
		signal d_r : sfixed(c_twiddle_int_width - 1 downto -c_twiddle_frac_width);
		-- stage B: products
		signal p_ac, p_bd, p_ad, p_bc : t_prod;
		-- stage C: full-precision components
		signal re_full, im_full : t_prod_sum;
		-- valid pipeline (A1, A2, B, C, D)
		signal v : std_logic_vector(1 to 5);
	begin

		PROC_PIPE: process(i_clk)
		begin
			if rising_edge(i_clk) then
				-- stage A1
				a_r1 <= i_sample.re;
				b_r1 <= i_sample.im;
				c_r1 <= i_twiddle.re;
				d_r1 <= i_twiddle.im;

				-- stage A2
				a_r <= a_r1;
				b_r <= b_r1;
				c_r <= c_r1;
				d_r <= d_r1;

				-- stage B: (a + jb)(c + jd) needs ac, bd, ad, bc
				p_ac <= a_r * c_r;
				p_bd <= b_r * d_r;
				p_ad <= a_r * d_r;
				p_bc <= b_r * c_r;

				-- stage C: full precision, no quantization yet
				re_full <= p_ac - p_bd;
				im_full <= p_ad + p_bc;

				-- stage D: the single wrap + half-up quantization (+half LSB
				-- then truncate), identical to the mr_fft_pkg "*" operator
				o_sample.re <= resize(re_full + c_half, o_sample.re, fixed_wrap, fixed_truncate);
				o_sample.im <= resize(im_full + c_half, o_sample.im, fixed_wrap, fixed_truncate);

				-- valid alongside
				if i_reset = '1' then
					v <= (others => '0');
				else
					v <= i_valid & v(1 to 4);
				end if;
			end if;
		end process PROC_PIPE;

		o_valid <= v(5);

	end generate GEN_PIPE;

end architecture rtl;
