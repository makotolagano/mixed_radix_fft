library ieee;
use ieee.std_logic_1164.all;
use ieee.fixed_pkg.all;
use ieee.fixed_float_types.all;

-- Custom packages
library work;
use work.mr_fft_pkg.all;

-- "18 bits in memory, wide in flight" (docs/datapath_width_convention.md):
-- inputs are the s2.16 memory word; the adder tree runs exact in the wide
-- word (s5.18 -- the s0 algorithmic shift keeps its bits), the coefficient
-- products are trimmed HALF-UP to the prod word (s6.22, +2^-23 then
-- truncate, foldable into the DSP post-adder), the t3 adds run exact in the
-- prod word, and each output applies the FIXED scaling shift derived from
-- the radix mode (radix-2 -> 1, radix-3 -> 2, radix-5 -> 3, from i_s0) fused
-- with the ONE half-up exit rounding back to the memory word (wrap).
-- Bit-exact vs MixedRadix_PreAdder_FXP(exit_round=True, internal_frac=22,
-- rounding='half_up', shift_bits=ceil(log2(radix))).
entity mr_fft_preadder is
	generic (
		G_CAPABILITY : natural := 2;
		-- false: combinational (original behavior). true: pipelined with
		-- preadder_latency(G_CAPABILITY, true) register stages; the datapath
		-- free-runs (no enables) -- the consumer qualifies outputs with its own
		-- valid pipeline of the same depth. Same arithmetic, bit-exact results.
		G_PIPELINE : boolean := false
	);
	port (
		i_clk : in  std_logic;

		i_x0 : in  t_cmplx;
		i_x1 : in  t_cmplx;
		i_x2 : in  t_cmplx;
		i_x3 : in  t_cmplx;
		i_x4 : in  t_cmplx;

		i_s0 : in  std_logic_vector(1 downto 0);
		i_s1 : in  std_logic;

		o_X0 : out t_cmplx;
		o_X1 : out t_cmplx;
		o_X2 : out t_cmplx;
		o_X3 : out t_cmplx;
		o_X4 : out t_cmplx
	);
end entity mr_fft_preadder;

architecture rtl of mr_fft_preadder is
	-- Type definitions
	type t_cmplx_wide_array is array (5 downto 0) of t_cmplx_wide;
	type t_cmplx_prod_array is array (5 downto 0) of t_cmplx_prod;

	-- raw coefficient product: wide (s5.18) x coeff (s2.16) real products
	-- combined with one growth bit; lives in the DSP M/P registers
	subtype t_raw is sfixed(c_fxp_int_wide_width + c_coeff_int_width
	                        downto -(c_fxp_frac_wide_width + c_coeff_frac_width));
	type t_cmplx_raw is record
		re : t_raw;
		im : t_raw;
	end record t_cmplx_raw;

	-- exit-round working word: one growth bit above the prod word for the
	-- output selection sums (bounded ~22.3, wrap never fires)
	subtype t_acc is sfixed(c_fxp_prod_int_width downto -c_fxp_prod_frac_width);

	-- half-up rounding constants (+half LSB then truncate):
	-- product trim to frac 22 -> +2^-23; exit round to frac 16 seen at the
	-- pre-shift scale -> +2^(shift-17)
	constant c_prod_half : sfixed(0 downto -(c_fxp_prod_frac_width + 1)) :=
		to_sfixed(2.0 ** (-(c_fxp_prod_frac_width + 1)), 0, -(c_fxp_prod_frac_width + 1));
	constant c_exit_half_1 : t_acc := to_sfixed(2.0 ** (1 - c_fxp_frac_width - 1),
	                                            c_fxp_prod_int_width, -c_fxp_prod_frac_width);
	constant c_exit_half_2 : t_acc := to_sfixed(2.0 ** (2 - c_fxp_frac_width - 1),
	                                            c_fxp_prod_int_width, -c_fxp_prod_frac_width);
	constant c_exit_half_3 : t_acc := to_sfixed(2.0 ** (3 - c_fxp_frac_width - 1),
	                                            c_fxp_prod_int_width, -c_fxp_prod_frac_width);

	signal input_x0_wide, input_x1_wide, input_x2_wide, input_x3_wide, input_x4_wide : t_cmplx_wide;

begin

	-- Resize input signals to wide fixed-point representation (exact)
	input_x0_wide <= resize(i_x0, input_x0_wide);
	input_x1_wide <= resize(i_x1, input_x1_wide);
	input_x2_wide <= resize(i_x2, input_x2_wide);
	input_x3_wide <= resize(i_x3, input_x3_wide);
	input_x4_wide <= resize(i_x4, input_x4_wide);

	GEN_235: if G_CAPABILITY = 2 and not G_PIPELINE generate
		PROC_CALC_235: process(input_x0_wide, input_x1_wide, input_x2_wide, input_x3_wide, input_x4_wide, i_s0, i_s1)
			variable t2_s0_mux_out                : t_cmplx_wide;
			variable t2_mul1_in                   : t_cmplx_coeff;
			variable r1, r2, r3, r4               : t_cmplx_raw;
			variable t0, t1                       : t_cmplx_wide_array;
			variable t2_add0, t2_add1             : t_cmplx_wide;
			variable t2, t3                       : t_cmplx_prod_array;
			variable v_shift                      : natural range 1 to 3;
			variable v_half                       : t_acc;
			variable v0, v1, v2, v3, v4           : t_acc;
		begin
			-- Stage 0 (exact, wide)
			t0(0) := input_x0_wide;
			t0(1) := input_x1_wide + input_x4_wide;
			t0(2) := input_x2_wide + input_x3_wide;
			t0(3) := input_x1_wide - input_x4_wide;
			t0(4) := input_x2_wide - input_x3_wide;

			-- Stage 1 (exact, wide)
			t1(0) := t0(0);
			t1(1) := t0(1) + t0(2);
			t1(2) := t0(1) - t0(2);
			t1(3) := t0(3);
			t1(4) := t0(4);
			t1(5) := t0(3) + t0(4);

			-- Stage 2: fabric sums (exact -- the s0 shift keeps its bits in
			-- the wide word's extra fraction bits) ...
			t2_add0 := t1(0) + t1(1);

			case i_s0 is
				when "00" => -- no shift
					t2_s0_mux_out := t1(1);
				when "01" => -- arithmetic shift right by 1
					t2_s0_mux_out := shift_right(t1(1), 1);
				when "10" => -- arithmetic shift right by 2
					t2_s0_mux_out := shift_right(t1(1), 2);
				when others =>
					t2_s0_mux_out := t1(1); -- default case, no shift
			end case;

			t2_add1 := t1(0) - t2_s0_mux_out;

			-- ... and the coefficient products, trimmed HALF-UP to the prod
			-- word (the +half constant folds into the DSP post-adder)
			case i_s1 is
				when '0' => -- mul with k6
					t2_mul1_in := c_k6;
				when '1' => -- mul with k2
					t2_mul1_in := c_k2;
				when others =>
					t2_mul1_in := c_k2; -- default case, no multiplication
			end case;

			r1.re := resize(t1(2).re * t2_mul1_in.re - t1(2).im * t2_mul1_in.im + c_prod_half, r1.re, fixed_wrap, fixed_truncate);
			r1.im := resize(t1(2).re * t2_mul1_in.im + t1(2).im * t2_mul1_in.re + c_prod_half, r1.im, fixed_wrap, fixed_truncate);
			r2.re := resize(t1(3).re * c_k3.re - t1(3).im * c_k3.im + c_prod_half, r2.re, fixed_wrap, fixed_truncate);
			r2.im := resize(t1(3).re * c_k3.im + t1(3).im * c_k3.re + c_prod_half, r2.im, fixed_wrap, fixed_truncate);
			r3.re := resize(t1(4).re * c_k5.re - t1(4).im * c_k5.im + c_prod_half, r3.re, fixed_wrap, fixed_truncate);
			r3.im := resize(t1(4).re * c_k5.im + t1(4).im * c_k5.re + c_prod_half, r3.im, fixed_wrap, fixed_truncate);
			r4.re := resize(t1(5).re * c_k4.re - t1(5).im * c_k4.im + c_prod_half, r4.re, fixed_wrap, fixed_truncate);
			r4.im := resize(t1(5).re * c_k4.im + t1(5).im * c_k4.re + c_prod_half, r4.im, fixed_wrap, fixed_truncate);

			-- trim (truncate the biased product) + j-mux
			case i_s1 is
				when '0' => -- mul with j
					t2(2).re := resize(-resize(r1.im, c_fxp_prod_int_width - 1, -c_fxp_prod_frac_width, fixed_wrap, fixed_truncate), t2(2).re, fixed_wrap, fixed_truncate);
					t2(2).im := resize(r1.re, t2(2).im, fixed_wrap, fixed_truncate);
				when '1' => -- passthrough
					t2(2).re := resize(r1.re, t2(2).re, fixed_wrap, fixed_truncate);
					t2(2).im := resize(r1.im, t2(2).im, fixed_wrap, fixed_truncate);
				when others =>
					t2(2).re := resize(r1.re, t2(2).re, fixed_wrap, fixed_truncate);
					t2(2).im := resize(r1.im, t2(2).im, fixed_wrap, fixed_truncate);
			end case;

			t2(3).re := resize(r2.re, t2(3).re, fixed_wrap, fixed_truncate);
			t2(3).im := resize(r2.im, t2(3).im, fixed_wrap, fixed_truncate);
			t2(4).re := resize(r3.re, t2(4).re, fixed_wrap, fixed_truncate);
			t2(4).im := resize(r3.im, t2(4).im, fixed_wrap, fixed_truncate);
			t2(5).re := resize(r4.re, t2(5).re, fixed_wrap, fixed_truncate);
			t2(5).im := resize(r4.im, t2(5).im, fixed_wrap, fixed_truncate);

			-- t2 sums carried into the prod word (exact resizes)
			t2(0).re := resize(t2_add0.re, t2(0).re, fixed_wrap, fixed_truncate);
			t2(0).im := resize(t2_add0.im, t2(0).im, fixed_wrap, fixed_truncate);
			t2(1).re := resize(t2_add1.re, t2(1).re, fixed_wrap, fixed_truncate);
			t2(1).im := resize(t2_add1.im, t2(1).im, fixed_wrap, fixed_truncate);

			-- Stage 3 (exact, prod word)
			t3(0) := t2(0);
			t3(5) := t2(1);
			t3(1).re := resize(t2(1).re + t2(2).re, t3(1).re, fixed_wrap, fixed_truncate);
			t3(1).im := resize(t2(1).im + t2(2).im, t3(1).im, fixed_wrap, fixed_truncate);
			t3(2).re := resize(t2(1).re - t2(2).re, t3(2).re, fixed_wrap, fixed_truncate);
			t3(2).im := resize(t2(1).im - t2(2).im, t3(2).im, fixed_wrap, fixed_truncate);
			t3(3).re := resize(t2(3).re + t2(5).re, t3(3).re, fixed_wrap, fixed_truncate);
			t3(3).im := resize(t2(3).im + t2(5).im, t3(3).im, fixed_wrap, fixed_truncate);
			t3(4).re := resize(t2(4).re + t2(5).re, t3(4).re, fixed_wrap, fixed_truncate);
			t3(4).im := resize(t2(4).im + t2(5).im, t3(4).im, fixed_wrap, fixed_truncate);

			-- Stage output: selection sums + fused scaling-shift/half-up exit
			-- rounding back to the memory word (the ONLY data rounding; wrap)
			case i_s0 is
				when "00"   => v_shift := 1; v_half := c_exit_half_1;  -- radix-2
				when "01"   => v_shift := 2; v_half := c_exit_half_2;  -- radix-3
				when others => v_shift := 3; v_half := c_exit_half_3;  -- radix-5
			end case;

			v0 := resize(t3(0).re, v0);
			case i_s0 is
				when "00"   => v1 := resize(t3(5).re, v1);
				when "01"   => v1 := resize(t3(1).re, v1);
				when others => v1 := resize(t3(1).re + t3(3).re, v1, fixed_wrap, fixed_truncate);
			end case;
			case i_s1 is
				when '0'    => v2 := resize(t3(2).re, v2);
				when others => v2 := resize(t3(2).re + t3(4).re, v2, fixed_wrap, fixed_truncate);
			end case;
			v3 := resize(t3(2).re - t3(4).re, v3, fixed_wrap, fixed_truncate);
			v4 := resize(t3(1).re - t3(3).re, v4, fixed_wrap, fixed_truncate);

			o_X0.re <= resize(shift_right(resize(v0 + v_half, v0, fixed_wrap, fixed_truncate), v_shift), o_X0.re, fixed_wrap, fixed_truncate);
			o_X1.re <= resize(shift_right(resize(v1 + v_half, v1, fixed_wrap, fixed_truncate), v_shift), o_X1.re, fixed_wrap, fixed_truncate);
			o_X2.re <= resize(shift_right(resize(v2 + v_half, v2, fixed_wrap, fixed_truncate), v_shift), o_X2.re, fixed_wrap, fixed_truncate);
			o_X3.re <= resize(shift_right(resize(v3 + v_half, v3, fixed_wrap, fixed_truncate), v_shift), o_X3.re, fixed_wrap, fixed_truncate);
			o_X4.re <= resize(shift_right(resize(v4 + v_half, v4, fixed_wrap, fixed_truncate), v_shift), o_X4.re, fixed_wrap, fixed_truncate);

			v0 := resize(t3(0).im, v0);
			case i_s0 is
				when "00"   => v1 := resize(t3(5).im, v1);
				when "01"   => v1 := resize(t3(1).im, v1);
				when others => v1 := resize(t3(1).im + t3(3).im, v1, fixed_wrap, fixed_truncate);
			end case;
			case i_s1 is
				when '0'    => v2 := resize(t3(2).im, v2);
				when others => v2 := resize(t3(2).im + t3(4).im, v2, fixed_wrap, fixed_truncate);
			end case;
			v3 := resize(t3(2).im - t3(4).im, v3, fixed_wrap, fixed_truncate);
			v4 := resize(t3(1).im - t3(3).im, v4, fixed_wrap, fixed_truncate);

			o_X0.im <= resize(shift_right(resize(v0 + v_half, v0, fixed_wrap, fixed_truncate), v_shift), o_X0.im, fixed_wrap, fixed_truncate);
			o_X1.im <= resize(shift_right(resize(v1 + v_half, v1, fixed_wrap, fixed_truncate), v_shift), o_X1.im, fixed_wrap, fixed_truncate);
			o_X2.im <= resize(shift_right(resize(v2 + v_half, v2, fixed_wrap, fixed_truncate), v_shift), o_X2.im, fixed_wrap, fixed_truncate);
			o_X3.im <= resize(shift_right(resize(v3 + v_half, v3, fixed_wrap, fixed_truncate), v_shift), o_X3.im, fixed_wrap, fixed_truncate);
			o_X4.im <= resize(shift_right(resize(v4 + v_half, v4, fixed_wrap, fixed_truncate), v_shift), o_X4.im, fixed_wrap, fixed_truncate);
		end process PROC_CALC_235;
	end generate GEN_235;

	GEN_23: if G_CAPABILITY = 1 and not G_PIPELINE generate
		PROC_CALC_23: process(input_x0_wide, input_x1_wide, input_x2_wide, i_s0)
			variable t1_s0_mux_out          : t_cmplx_wide;
			variable r1                     : t_cmplx_raw;
			variable t0                     : t_cmplx_wide_array;
			variable t1_add0, t1_add1       : t_cmplx_wide;
			variable t1_2                   : t_cmplx_prod;
			variable t2                     : t_cmplx_prod_array;
			variable v_shift                : natural range 1 to 2;
			variable v_half                 : t_acc;
			variable v0, v1, v2             : t_acc;
		begin
			-- Stage 0 (exact, wide)
			t0(0) := input_x0_wide;
			t0(1) := input_x1_wide + input_x2_wide;
			t0(2) := input_x1_wide - input_x2_wide;

			-- Stage 1: fabric sums (exact) + trimmed k6 product
			t1_add0 := t0(0) + t0(1);

			case i_s0(0) is
				when '0' => -- no shift
					t1_s0_mux_out := t0(1);
				when '1' => -- arithmetic shift right by 1
					t1_s0_mux_out := shift_right(t0(1), 1);
				when others =>
					t1_s0_mux_out := t0(1); -- default case, no shift
			end case;

			t1_add1 := t0(0) - t1_s0_mux_out;

			r1.re := resize(t0(2).re * c_k6.re - t0(2).im * c_k6.im + c_prod_half, r1.re, fixed_wrap, fixed_truncate);
			r1.im := resize(t0(2).re * c_k6.im + t0(2).im * c_k6.re + c_prod_half, r1.im, fixed_wrap, fixed_truncate);

			-- trim + multiplication by j
			t1_2.re := resize(-resize(r1.im, c_fxp_prod_int_width - 1, -c_fxp_prod_frac_width, fixed_wrap, fixed_truncate), t1_2.re, fixed_wrap, fixed_truncate);
			t1_2.im := resize(r1.re, t1_2.im, fixed_wrap, fixed_truncate);

			-- Stage 2 (exact, prod word)
			t2(0).re := resize(t1_add0.re, t2(0).re, fixed_wrap, fixed_truncate);
			t2(0).im := resize(t1_add0.im, t2(0).im, fixed_wrap, fixed_truncate);
			t2(3).re := resize(t1_add1.re, t2(3).re, fixed_wrap, fixed_truncate);
			t2(3).im := resize(t1_add1.im, t2(3).im, fixed_wrap, fixed_truncate);
			t2(1).re := resize(t2(3).re + t1_2.re, t2(1).re, fixed_wrap, fixed_truncate);
			t2(1).im := resize(t2(3).im + t1_2.im, t2(1).im, fixed_wrap, fixed_truncate);
			t2(2).re := resize(t2(3).re - t1_2.re, t2(2).re, fixed_wrap, fixed_truncate);
			t2(2).im := resize(t2(3).im - t1_2.im, t2(2).im, fixed_wrap, fixed_truncate);

			-- Stage output: fused scaling-shift/half-up exit rounding
			case i_s0(0) is
				when '0'    => v_shift := 1; v_half := c_exit_half_1;  -- radix-2
				when others => v_shift := 2; v_half := c_exit_half_2;  -- radix-3
			end case;

			v0 := resize(t2(0).re, v0);
			case i_s0(0) is
				when '0'    => v1 := resize(t2(3).re, v1);
				when others => v1 := resize(t2(1).re, v1);
			end case;
			v2 := resize(t2(2).re, v2);
			o_X0.re <= resize(shift_right(resize(v0 + v_half, v0, fixed_wrap, fixed_truncate), v_shift), o_X0.re, fixed_wrap, fixed_truncate);
			o_X1.re <= resize(shift_right(resize(v1 + v_half, v1, fixed_wrap, fixed_truncate), v_shift), o_X1.re, fixed_wrap, fixed_truncate);
			o_X2.re <= resize(shift_right(resize(v2 + v_half, v2, fixed_wrap, fixed_truncate), v_shift), o_X2.re, fixed_wrap, fixed_truncate);

			v0 := resize(t2(0).im, v0);
			case i_s0(0) is
				when '0'    => v1 := resize(t2(3).im, v1);
				when others => v1 := resize(t2(1).im, v1);
			end case;
			v2 := resize(t2(2).im, v2);
			o_X0.im <= resize(shift_right(resize(v0 + v_half, v0, fixed_wrap, fixed_truncate), v_shift), o_X0.im, fixed_wrap, fixed_truncate);
			o_X1.im <= resize(shift_right(resize(v1 + v_half, v1, fixed_wrap, fixed_truncate), v_shift), o_X1.im, fixed_wrap, fixed_truncate);
			o_X2.im <= resize(shift_right(resize(v2 + v_half, v2, fixed_wrap, fixed_truncate), v_shift), o_X2.im, fixed_wrap, fixed_truncate);
		end process PROC_CALC_23;
	end generate GEN_23;

	GEN_2: if G_CAPABILITY = 0 and not G_PIPELINE generate
		PROC_CALC_2: process(input_x0_wide, input_x1_wide)
			variable t1                 : t_cmplx_wide_array;
			variable v0, v1             : t_acc;
		begin
			-- Stage 1 (exact, wide)
			t1(0) := input_x0_wide + input_x1_wide;
			t1(1) := input_x0_wide - input_x1_wide;

			-- Stage output: radix-2 scaling shift (1) fused with the half-up
			-- exit rounding
			v0 := resize(t1(0).re, v0);
			v1 := resize(t1(1).re, v1);
			o_X0.re <= resize(shift_right(resize(v0 + c_exit_half_1, v0, fixed_wrap, fixed_truncate), 1), o_X0.re, fixed_wrap, fixed_truncate);
			o_X1.re <= resize(shift_right(resize(v1 + c_exit_half_1, v1, fixed_wrap, fixed_truncate), 1), o_X1.re, fixed_wrap, fixed_truncate);

			v0 := resize(t1(0).im, v0);
			v1 := resize(t1(1).im, v1);
			o_X0.im <= resize(shift_right(resize(v0 + c_exit_half_1, v0, fixed_wrap, fixed_truncate), 1), o_X0.im, fixed_wrap, fixed_truncate);
			o_X1.im <= resize(shift_right(resize(v1 + c_exit_half_1, v1, fixed_wrap, fixed_truncate), 1), o_X1.im, fixed_wrap, fixed_truncate);
		end process PROC_CALC_2;
	end generate GEN_2;

	-- ------------------------------------------------------------------
	-- Pipelined variants: identical arithmetic to the combinational ones,
	-- with register cuts per preadder_latency(). Registers free-run.
	-- The biased raw products (product + half-LSB, i.e. the trim rounding
	-- constant on the DSP post-adder) are registered directly so Vivado
	-- absorbs the registers as DSP M/P registers.
	-- ------------------------------------------------------------------

	GEN_235_PIPE: if G_CAPABILITY = 2 and G_PIPELINE generate  -- latency 3
		signal t1_r : t_cmplx_wide_array;
		signal t2_add0_r, t2_add1_r : t_cmplx_wide;
		signal r1_r, r2_r, r3_r, r4_r : t_cmplx_raw;
		signal t3_r : t_cmplx_prod_array;

		attribute use_dsp : string;
		attribute use_dsp of r1_r, r2_r, r3_r, r4_r : signal is "yes";
	begin
		-- cut 1: after the two input adder levels (t0, t1)
		PROC_STAGE_A: process(i_clk)
			variable t0, t1 : t_cmplx_wide_array;
		begin
			if rising_edge(i_clk) then
				t0(0) := input_x0_wide;
				t0(1) := input_x1_wide + input_x4_wide;
				t0(2) := input_x2_wide + input_x3_wide;
				t0(3) := input_x1_wide - input_x4_wide;
				t0(4) := input_x2_wide - input_x3_wide;

				t1(0) := t0(0);
				t1(1) := t0(1) + t0(2);
				t1(2) := t0(1) - t0(2);
				t1(3) := t0(3);
				t1(4) := t0(4);
				t1(5) := t0(3) + t0(4);

				t1_r <= t1;
			end if;
		end process PROC_STAGE_A;

		-- cut 2: the t2 sums plus the BIASED raw products (product + trim
		-- half-LSB on the DSP post-adder -> registers absorb as M/P regs)
		PROC_STAGE_B: process(i_clk)
			variable t2_s0_mux_out : t_cmplx_wide;
			variable t2_mul1_in : t_cmplx_coeff;
		begin
			if rising_edge(i_clk) then
				t2_add0_r <= t1_r(0) + t1_r(1);

				case i_s0 is
					when "00" => -- no shift
						t2_s0_mux_out := t1_r(1);
					when "01" => -- arithmetic shift right by 1
						t2_s0_mux_out := shift_right(t1_r(1), 1);
					when "10" => -- arithmetic shift right by 2
						t2_s0_mux_out := shift_right(t1_r(1), 2);
					when others =>
						t2_s0_mux_out := t1_r(1); -- default case, no shift
				end case;

				t2_add1_r <= t1_r(0) - t2_s0_mux_out;

				case i_s1 is
					when '0' => -- mul with k6
						t2_mul1_in := c_k6;
					when '1' => -- mul with k2
						t2_mul1_in := c_k2;
					when others =>
						t2_mul1_in := c_k2; -- default case, no multiplication
				end case;

				r1_r.re <= resize(t1_r(2).re * t2_mul1_in.re - t1_r(2).im * t2_mul1_in.im + c_prod_half, r1_r.re, fixed_wrap, fixed_truncate);
				r1_r.im <= resize(t1_r(2).re * t2_mul1_in.im + t1_r(2).im * t2_mul1_in.re + c_prod_half, r1_r.im, fixed_wrap, fixed_truncate);
				r2_r.re <= resize(t1_r(3).re * c_k3.re - t1_r(3).im * c_k3.im + c_prod_half, r2_r.re, fixed_wrap, fixed_truncate);
				r2_r.im <= resize(t1_r(3).re * c_k3.im + t1_r(3).im * c_k3.re + c_prod_half, r2_r.im, fixed_wrap, fixed_truncate);
				r3_r.re <= resize(t1_r(4).re * c_k5.re - t1_r(4).im * c_k5.im + c_prod_half, r3_r.re, fixed_wrap, fixed_truncate);
				r3_r.im <= resize(t1_r(4).re * c_k5.im + t1_r(4).im * c_k5.re + c_prod_half, r3_r.im, fixed_wrap, fixed_truncate);
				r4_r.re <= resize(t1_r(5).re * c_k4.re - t1_r(5).im * c_k4.im + c_prod_half, r4_r.re, fixed_wrap, fixed_truncate);
				r4_r.im <= resize(t1_r(5).re * c_k4.im + t1_r(5).im * c_k4.re + c_prod_half, r4_r.im, fixed_wrap, fixed_truncate);
			end if;
		end process PROC_STAGE_B;

		-- cut 3: trim slice + j-mux + t3 adds (prod word)
		PROC_STAGE_C: process(i_clk)
			variable t2, t3 : t_cmplx_prod_array;
		begin
			if rising_edge(i_clk) then
				case i_s1 is
					when '0' => -- mul with j
						t2(2).re := resize(-resize(r1_r.im, c_fxp_prod_int_width - 1, -c_fxp_prod_frac_width, fixed_wrap, fixed_truncate), t2(2).re, fixed_wrap, fixed_truncate);
						t2(2).im := resize(r1_r.re, t2(2).im, fixed_wrap, fixed_truncate);
					when '1' => -- passthrough
						t2(2).re := resize(r1_r.re, t2(2).re, fixed_wrap, fixed_truncate);
						t2(2).im := resize(r1_r.im, t2(2).im, fixed_wrap, fixed_truncate);
					when others =>
						t2(2).re := resize(r1_r.re, t2(2).re, fixed_wrap, fixed_truncate);
						t2(2).im := resize(r1_r.im, t2(2).im, fixed_wrap, fixed_truncate);
				end case;

				t2(3).re := resize(r2_r.re, t2(3).re, fixed_wrap, fixed_truncate);
				t2(3).im := resize(r2_r.im, t2(3).im, fixed_wrap, fixed_truncate);
				t2(4).re := resize(r3_r.re, t2(4).re, fixed_wrap, fixed_truncate);
				t2(4).im := resize(r3_r.im, t2(4).im, fixed_wrap, fixed_truncate);
				t2(5).re := resize(r4_r.re, t2(5).re, fixed_wrap, fixed_truncate);
				t2(5).im := resize(r4_r.im, t2(5).im, fixed_wrap, fixed_truncate);

				t2(0).re := resize(t2_add0_r.re, t2(0).re, fixed_wrap, fixed_truncate);
				t2(0).im := resize(t2_add0_r.im, t2(0).im, fixed_wrap, fixed_truncate);
				t2(1).re := resize(t2_add1_r.re, t2(1).re, fixed_wrap, fixed_truncate);
				t2(1).im := resize(t2_add1_r.im, t2(1).im, fixed_wrap, fixed_truncate);

				t3(0) := t2(0);
				t3(5) := t2(1);
				t3(1).re := resize(t2(1).re + t2(2).re, t3(1).re, fixed_wrap, fixed_truncate);
				t3(1).im := resize(t2(1).im + t2(2).im, t3(1).im, fixed_wrap, fixed_truncate);
				t3(2).re := resize(t2(1).re - t2(2).re, t3(2).re, fixed_wrap, fixed_truncate);
				t3(2).im := resize(t2(1).im - t2(2).im, t3(2).im, fixed_wrap, fixed_truncate);
				t3(3).re := resize(t2(3).re + t2(5).re, t3(3).re, fixed_wrap, fixed_truncate);
				t3(3).im := resize(t2(3).im + t2(5).im, t3(3).im, fixed_wrap, fixed_truncate);
				t3(4).re := resize(t2(4).re + t2(5).re, t3(4).re, fixed_wrap, fixed_truncate);
				t3(4).im := resize(t2(4).im + t2(5).im, t3(4).im, fixed_wrap, fixed_truncate);

				t3_r <= t3;
			end if;
		end process PROC_STAGE_C;

		-- output stage (combinational): selection sums + fused shift/round
		PROC_STAGE_D: process(t3_r, i_s0, i_s1)
			variable v_shift            : natural range 1 to 3;
			variable v_half             : t_acc;
			variable v0, v1, v2, v3, v4 : t_acc;
		begin
			case i_s0 is
				when "00"   => v_shift := 1; v_half := c_exit_half_1;  -- radix-2
				when "01"   => v_shift := 2; v_half := c_exit_half_2;  -- radix-3
				when others => v_shift := 3; v_half := c_exit_half_3;  -- radix-5
			end case;

			v0 := resize(t3_r(0).re, v0);
			case i_s0 is
				when "00"   => v1 := resize(t3_r(5).re, v1);
				when "01"   => v1 := resize(t3_r(1).re, v1);
				when others => v1 := resize(t3_r(1).re + t3_r(3).re, v1, fixed_wrap, fixed_truncate);
			end case;
			case i_s1 is
				when '0'    => v2 := resize(t3_r(2).re, v2);
				when others => v2 := resize(t3_r(2).re + t3_r(4).re, v2, fixed_wrap, fixed_truncate);
			end case;
			v3 := resize(t3_r(2).re - t3_r(4).re, v3, fixed_wrap, fixed_truncate);
			v4 := resize(t3_r(1).re - t3_r(3).re, v4, fixed_wrap, fixed_truncate);

			o_X0.re <= resize(shift_right(resize(v0 + v_half, v0, fixed_wrap, fixed_truncate), v_shift), o_X0.re, fixed_wrap, fixed_truncate);
			o_X1.re <= resize(shift_right(resize(v1 + v_half, v1, fixed_wrap, fixed_truncate), v_shift), o_X1.re, fixed_wrap, fixed_truncate);
			o_X2.re <= resize(shift_right(resize(v2 + v_half, v2, fixed_wrap, fixed_truncate), v_shift), o_X2.re, fixed_wrap, fixed_truncate);
			o_X3.re <= resize(shift_right(resize(v3 + v_half, v3, fixed_wrap, fixed_truncate), v_shift), o_X3.re, fixed_wrap, fixed_truncate);
			o_X4.re <= resize(shift_right(resize(v4 + v_half, v4, fixed_wrap, fixed_truncate), v_shift), o_X4.re, fixed_wrap, fixed_truncate);

			v0 := resize(t3_r(0).im, v0);
			case i_s0 is
				when "00"   => v1 := resize(t3_r(5).im, v1);
				when "01"   => v1 := resize(t3_r(1).im, v1);
				when others => v1 := resize(t3_r(1).im + t3_r(3).im, v1, fixed_wrap, fixed_truncate);
			end case;
			case i_s1 is
				when '0'    => v2 := resize(t3_r(2).im, v2);
				when others => v2 := resize(t3_r(2).im + t3_r(4).im, v2, fixed_wrap, fixed_truncate);
			end case;
			v3 := resize(t3_r(2).im - t3_r(4).im, v3, fixed_wrap, fixed_truncate);
			v4 := resize(t3_r(1).im - t3_r(3).im, v4, fixed_wrap, fixed_truncate);

			o_X0.im <= resize(shift_right(resize(v0 + v_half, v0, fixed_wrap, fixed_truncate), v_shift), o_X0.im, fixed_wrap, fixed_truncate);
			o_X1.im <= resize(shift_right(resize(v1 + v_half, v1, fixed_wrap, fixed_truncate), v_shift), o_X1.im, fixed_wrap, fixed_truncate);
			o_X2.im <= resize(shift_right(resize(v2 + v_half, v2, fixed_wrap, fixed_truncate), v_shift), o_X2.im, fixed_wrap, fixed_truncate);
			o_X3.im <= resize(shift_right(resize(v3 + v_half, v3, fixed_wrap, fixed_truncate), v_shift), o_X3.im, fixed_wrap, fixed_truncate);
			o_X4.im <= resize(shift_right(resize(v4 + v_half, v4, fixed_wrap, fixed_truncate), v_shift), o_X4.im, fixed_wrap, fixed_truncate);
		end process PROC_STAGE_D;
	end generate GEN_235_PIPE;

	GEN_23_PIPE: if G_CAPABILITY = 1 and G_PIPELINE generate  -- latency 3
		signal t0_r : t_cmplx_wide_array;
		signal t1_add0_r, t1_add1_r : t_cmplx_wide;
		signal r1_r : t_cmplx_raw;
		signal t2_r : t_cmplx_prod_array;

		attribute use_dsp : string;
		attribute use_dsp of r1_r : signal is "yes";
	begin
		-- cut 1: after the input adder level (t0)
		PROC_STAGE_A: process(i_clk)
			variable t0 : t_cmplx_wide_array;
		begin
			if rising_edge(i_clk) then
				t0(0) := input_x0_wide;
				t0(1) := input_x1_wide + input_x2_wide;
				t0(2) := input_x1_wide - input_x2_wide;
				t0_r <= t0;
			end if;
		end process PROC_STAGE_A;

		-- cut 2: the t1 sums plus the BIASED raw k6 product
		PROC_STAGE_B: process(i_clk)
			variable t1_s0_mux_out : t_cmplx_wide;
		begin
			if rising_edge(i_clk) then
				t1_add0_r <= t0_r(0) + t0_r(1);

				case i_s0(0) is
					when '0' => -- no shift
						t1_s0_mux_out := t0_r(1);
					when '1' => -- arithmetic shift right by 1
						t1_s0_mux_out := shift_right(t0_r(1), 1);
					when others =>
						t1_s0_mux_out := t0_r(1); -- default case, no shift
				end case;

				t1_add1_r <= t0_r(0) - t1_s0_mux_out;

				r1_r.re <= resize(t0_r(2).re * c_k6.re - t0_r(2).im * c_k6.im + c_prod_half, r1_r.re, fixed_wrap, fixed_truncate);
				r1_r.im <= resize(t0_r(2).re * c_k6.im + t0_r(2).im * c_k6.re + c_prod_half, r1_r.im, fixed_wrap, fixed_truncate);
			end if;
		end process PROC_STAGE_B;

		-- cut 3: trim slice + multiplication by j + t2 adds (prod word)
		PROC_STAGE_C: process(i_clk)
			variable t1_2 : t_cmplx_prod;
			variable t2 : t_cmplx_prod_array;
		begin
			if rising_edge(i_clk) then
				t1_2.re := resize(-resize(r1_r.im, c_fxp_prod_int_width - 1, -c_fxp_prod_frac_width, fixed_wrap, fixed_truncate), t1_2.re, fixed_wrap, fixed_truncate);
				t1_2.im := resize(r1_r.re, t1_2.im, fixed_wrap, fixed_truncate);

				t2(0).re := resize(t1_add0_r.re, t2(0).re, fixed_wrap, fixed_truncate);
				t2(0).im := resize(t1_add0_r.im, t2(0).im, fixed_wrap, fixed_truncate);
				t2(3).re := resize(t1_add1_r.re, t2(3).re, fixed_wrap, fixed_truncate);
				t2(3).im := resize(t1_add1_r.im, t2(3).im, fixed_wrap, fixed_truncate);
				t2(1).re := resize(t2(3).re + t1_2.re, t2(1).re, fixed_wrap, fixed_truncate);
				t2(1).im := resize(t2(3).im + t1_2.im, t2(1).im, fixed_wrap, fixed_truncate);
				t2(2).re := resize(t2(3).re - t1_2.re, t2(2).re, fixed_wrap, fixed_truncate);
				t2(2).im := resize(t2(3).im - t1_2.im, t2(2).im, fixed_wrap, fixed_truncate);

				t2_r <= t2;
			end if;
		end process PROC_STAGE_C;

		-- output stage (combinational): selection + fused shift/round
		PROC_STAGE_D: process(t2_r, i_s0)
			variable v_shift    : natural range 1 to 2;
			variable v_half     : t_acc;
			variable v0, v1, v2 : t_acc;
		begin
			case i_s0(0) is
				when '0'    => v_shift := 1; v_half := c_exit_half_1;  -- radix-2
				when others => v_shift := 2; v_half := c_exit_half_2;  -- radix-3
			end case;

			v0 := resize(t2_r(0).re, v0);
			case i_s0(0) is
				when '0'    => v1 := resize(t2_r(3).re, v1);
				when others => v1 := resize(t2_r(1).re, v1);
			end case;
			v2 := resize(t2_r(2).re, v2);
			o_X0.re <= resize(shift_right(resize(v0 + v_half, v0, fixed_wrap, fixed_truncate), v_shift), o_X0.re, fixed_wrap, fixed_truncate);
			o_X1.re <= resize(shift_right(resize(v1 + v_half, v1, fixed_wrap, fixed_truncate), v_shift), o_X1.re, fixed_wrap, fixed_truncate);
			o_X2.re <= resize(shift_right(resize(v2 + v_half, v2, fixed_wrap, fixed_truncate), v_shift), o_X2.re, fixed_wrap, fixed_truncate);

			v0 := resize(t2_r(0).im, v0);
			case i_s0(0) is
				when '0'    => v1 := resize(t2_r(3).im, v1);
				when others => v1 := resize(t2_r(1).im, v1);
			end case;
			v2 := resize(t2_r(2).im, v2);
			o_X0.im <= resize(shift_right(resize(v0 + v_half, v0, fixed_wrap, fixed_truncate), v_shift), o_X0.im, fixed_wrap, fixed_truncate);
			o_X1.im <= resize(shift_right(resize(v1 + v_half, v1, fixed_wrap, fixed_truncate), v_shift), o_X1.im, fixed_wrap, fixed_truncate);
			o_X2.im <= resize(shift_right(resize(v2 + v_half, v2, fixed_wrap, fixed_truncate), v_shift), o_X2.im, fixed_wrap, fixed_truncate);
		end process PROC_STAGE_D;
	end generate GEN_23_PIPE;

	GEN_2_PIPE: if G_CAPABILITY = 0 and G_PIPELINE generate  -- latency 1
		PROC_STAGE_A: process(i_clk)
			variable t1     : t_cmplx_wide_array;
			variable v0, v1 : t_acc;
		begin
			if rising_edge(i_clk) then
				t1(0) := input_x0_wide + input_x1_wide;
				t1(1) := input_x0_wide - input_x1_wide;

				v0 := resize(t1(0).re, v0);
				v1 := resize(t1(1).re, v1);
				o_X0.re <= resize(shift_right(resize(v0 + c_exit_half_1, v0, fixed_wrap, fixed_truncate), 1), o_X0.re, fixed_wrap, fixed_truncate);
				o_X1.re <= resize(shift_right(resize(v1 + c_exit_half_1, v1, fixed_wrap, fixed_truncate), 1), o_X1.re, fixed_wrap, fixed_truncate);

				v0 := resize(t1(0).im, v0);
				v1 := resize(t1(1).im, v1);
				o_X0.im <= resize(shift_right(resize(v0 + c_exit_half_1, v0, fixed_wrap, fixed_truncate), 1), o_X0.im, fixed_wrap, fixed_truncate);
				o_X1.im <= resize(shift_right(resize(v1 + c_exit_half_1, v1, fixed_wrap, fixed_truncate), 1), o_X1.im, fixed_wrap, fixed_truncate);
			end if;
		end process PROC_STAGE_A;
	end generate GEN_2_PIPE;

end architecture rtl;
