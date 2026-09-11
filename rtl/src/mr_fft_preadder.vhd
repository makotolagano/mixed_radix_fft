library ieee;
use ieee.std_logic_1164.all;
use ieee.fixed_pkg.all;
use ieee.fixed_float_types.all;

-- Custom packages
library work;
use work.mr_fft_pkg.all;

-- Reconfigurable preadder (butterfly) for radix 2/3/5.
-- inputs are the s2.16 word, the adder tree runs in a wider word, coefficient
-- products are trimmed half-up, and every output gets the radix scaling shift
-- (1/2/3 bits) fused with one half-up rounding back to s2.16.
entity mr_fft_preadder is
	generic (
		G_CAPABILITY : natural := 2;
		-- true: register cuts per preadder_latency, no enables. the consumer runs its
		-- own valid pipeline. same arithmetic either way.
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

	-- one output component (t_cmplx.re/.im)
	subtype t_comp is sfixed(c_fxp_int_width - 1 downto -c_fxp_frac_width);

	-- raw coefficient product, lives in the DSP M/P registers
	subtype t_raw is sfixed(c_fxp_int_wide_width + c_coeff_int_width
	                        downto -(c_fxp_frac_wide_width + c_coeff_frac_width));
	type t_cmplx_raw is record
		re : t_raw;
		im : t_raw;
	end record t_cmplx_raw;

	-- working word for the output sums, one growth bit
	subtype t_acc is sfixed(c_fxp_prod_int_width downto -c_fxp_prod_frac_width);

	-- half-up rounding constants (+half LSB then truncate)
	constant c_prod_half : sfixed(0 downto -(c_fxp_prod_frac_width + 1)) :=
		to_sfixed(2.0 ** (-(c_fxp_prod_frac_width + 1)), 0, -(c_fxp_prod_frac_width + 1));
	constant c_exit_half_1 : t_acc := to_sfixed(2.0 ** (1 - c_fxp_frac_width - 1),
	                                            c_fxp_prod_int_width, -c_fxp_prod_frac_width);
	constant c_exit_half_2 : t_acc := to_sfixed(2.0 ** (2 - c_fxp_frac_width - 1),
	                                            c_fxp_prod_int_width, -c_fxp_prod_frac_width);
	constant c_exit_half_3 : t_acc := to_sfixed(2.0 ** (3 - c_fxp_frac_width - 1),
	                                            c_fxp_prod_int_width, -c_fxp_prod_frac_width);

	-- scaling shift plus half-up rounding back to s2.16. shift amounts are literals,
	-- xsim 2022.2 crashes on shift_right with a variable amount.
	-- s0: "00" radix-2 (shift 1), "01" radix-3 (2), else radix-5 (3)
	function f_exit_round(value : t_acc; s0 : std_logic_vector(1 downto 0)) return t_comp is
		variable v : t_acc;
	begin
		case s0 is
			when "00" =>
				v := resize(value + c_exit_half_1, v, fixed_wrap, fixed_truncate);
				v := shift_right(v, 1);
			when "01" =>
				v := resize(value + c_exit_half_2, v, fixed_wrap, fixed_truncate);
				v := shift_right(v, 2);
			when others =>
				v := resize(value + c_exit_half_3, v, fixed_wrap, fixed_truncate);
				v := shift_right(v, 3);
		end case;
		return resize(v, c_fxp_int_width - 1, -c_fxp_frac_width, fixed_wrap, fixed_truncate);
	end function f_exit_round;

	signal input_x0_wide, input_x1_wide, input_x2_wide, input_x3_wide, input_x4_wide : t_cmplx_wide;

begin

	-- widen the inputs (exact)
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
			variable v0, v1, v2, v3, v4           : t_acc;
		begin
			-- stage 0
			t0(0) := input_x0_wide;
			t0(1) := input_x1_wide + input_x4_wide;
			t0(2) := input_x2_wide + input_x3_wide;
			t0(3) := input_x1_wide - input_x4_wide;
			t0(4) := input_x2_wide - input_x3_wide;

			-- stage 1
			t1(0) := t0(0);
			t1(1) := t0(1) + t0(2);
			t1(2) := t0(1) - t0(2);
			t1(3) := t0(3);
			t1(4) := t0(4);
			t1(5) := t0(3) + t0(4);

			-- stage 2: sums (the s0 shift keeps its bits in the wide word)
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

			-- coefficient products, trimmed half-up to the prod word
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

			-- trim and multiply by j
			case i_s1 is
				when '0' => -- mul with j
					t2(2).im := resize(r1.re, t2(2).im, fixed_wrap, fixed_truncate);
					t2(2).re := resize(r1.im, t2(2).re, fixed_wrap, fixed_truncate);
					t2(2).re := resize(-t2(2).re, t2(2).re, fixed_wrap, fixed_truncate);
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

			-- carry the sums into the prod word
			t2(0).re := resize(t2_add0.re, t2(0).re, fixed_wrap, fixed_truncate);
			t2(0).im := resize(t2_add0.im, t2(0).im, fixed_wrap, fixed_truncate);
			t2(1).re := resize(t2_add1.re, t2(1).re, fixed_wrap, fixed_truncate);
			t2(1).im := resize(t2_add1.im, t2(1).im, fixed_wrap, fixed_truncate);

			-- stage 3
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

			-- output: selection sums, then shift and round
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

			o_X0.re <= f_exit_round(v0, i_s0);
			o_X1.re <= f_exit_round(v1, i_s0);
			o_X2.re <= f_exit_round(v2, i_s0);
			o_X3.re <= f_exit_round(v3, i_s0);
			o_X4.re <= f_exit_round(v4, i_s0);

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

			o_X0.im <= f_exit_round(v0, i_s0);
			o_X1.im <= f_exit_round(v1, i_s0);
			o_X2.im <= f_exit_round(v2, i_s0);
			o_X3.im <= f_exit_round(v3, i_s0);
			o_X4.im <= f_exit_round(v4, i_s0);
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
			variable v0, v1, v2             : t_acc;
		begin
			-- stage 0
			t0(0) := input_x0_wide;
			t0(1) := input_x1_wide + input_x2_wide;
			t0(2) := input_x1_wide - input_x2_wide;

			-- stage 1: sums and the k6 product
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

			-- trim and multiply by j
			t1_2.im := resize(r1.re, t1_2.im, fixed_wrap, fixed_truncate);
			t1_2.re := resize(r1.im, t1_2.re, fixed_wrap, fixed_truncate);
			t1_2.re := resize(-t1_2.re, t1_2.re, fixed_wrap, fixed_truncate);

			-- stage 2
			t2(0).re := resize(t1_add0.re, t2(0).re, fixed_wrap, fixed_truncate);
			t2(0).im := resize(t1_add0.im, t2(0).im, fixed_wrap, fixed_truncate);
			t2(3).re := resize(t1_add1.re, t2(3).re, fixed_wrap, fixed_truncate);
			t2(3).im := resize(t1_add1.im, t2(3).im, fixed_wrap, fixed_truncate);
			t2(1).re := resize(t2(3).re + t1_2.re, t2(1).re, fixed_wrap, fixed_truncate);
			t2(1).im := resize(t2(3).im + t1_2.im, t2(1).im, fixed_wrap, fixed_truncate);
			t2(2).re := resize(t2(3).re - t1_2.re, t2(2).re, fixed_wrap, fixed_truncate);
			t2(2).im := resize(t2(3).im - t1_2.im, t2(2).im, fixed_wrap, fixed_truncate);

			-- output: selection, shift and round
			v0 := resize(t2(0).re, v0);
			case i_s0(0) is
				when '0'    => v1 := resize(t2(3).re, v1);
				when others => v1 := resize(t2(1).re, v1);
			end case;
			v2 := resize(t2(2).re, v2);
			o_X0.re <= f_exit_round(v0, '0' & i_s0(0));
			o_X1.re <= f_exit_round(v1, '0' & i_s0(0));
			o_X2.re <= f_exit_round(v2, '0' & i_s0(0));

			v0 := resize(t2(0).im, v0);
			case i_s0(0) is
				when '0'    => v1 := resize(t2(3).im, v1);
				when others => v1 := resize(t2(1).im, v1);
			end case;
			v2 := resize(t2(2).im, v2);
			o_X0.im <= f_exit_round(v0, '0' & i_s0(0));
			o_X1.im <= f_exit_round(v1, '0' & i_s0(0));
			o_X2.im <= f_exit_round(v2, '0' & i_s0(0));
		end process PROC_CALC_23;
	end generate GEN_23;

	GEN_2: if G_CAPABILITY = 0 and not G_PIPELINE generate
		PROC_CALC_2: process(input_x0_wide, input_x1_wide)
			variable t1     : t_cmplx_wide_array;
			variable v0, v1 : t_acc;
		begin
			-- stage 1
			t1(0) := input_x0_wide + input_x1_wide;
			t1(1) := input_x0_wide - input_x1_wide;

			-- output: shift by 1 and round
			v0 := resize(t1(0).re, v0);
			v1 := resize(t1(1).re, v1);
			o_X0.re <= f_exit_round(v0, "00");
			o_X1.re <= f_exit_round(v1, "00");

			v0 := resize(t1(0).im, v0);
			v1 := resize(t1(1).im, v1);
			o_X0.im <= f_exit_round(v0, "00");
			o_X1.im <= f_exit_round(v1, "00");
		end process PROC_CALC_2;
	end generate GEN_2;

	-- pipelined variants: same arithmetic, register cuts per preadder_latency.
	-- the raw products (with the rounding bias) are registered directly so
	-- Vivado absorbs them into the DSP M/P registers.

	GEN_235_PIPE: if G_CAPABILITY = 2 and G_PIPELINE generate  -- latency 3
		signal t1_r : t_cmplx_wide_array;
		signal t2_add0_r, t2_add1_r : t_cmplx_wide;
		signal r1_r, r2_r, r3_r, r4_r : t_cmplx_raw;
		signal t3_r : t_cmplx_prod_array;

		attribute use_dsp : string;
		attribute use_dsp of r1_r, r2_r, r3_r, r4_r : signal is "yes";
	begin
		-- cut 1: after the two input adder levels
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

		-- cut 2: t2 sums and the biased raw products
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

		-- cut 3: trim, j-mux and t3 adds
		PROC_STAGE_C: process(i_clk)
			variable t2, t3 : t_cmplx_prod_array;
		begin
			if rising_edge(i_clk) then
				case i_s1 is
					when '0' => -- mul with j
						t2(2).im := resize(r1_r.re, t2(2).im, fixed_wrap, fixed_truncate);
						t2(2).re := resize(r1_r.im, t2(2).re, fixed_wrap, fixed_truncate);
						t2(2).re := resize(-t2(2).re, t2(2).re, fixed_wrap, fixed_truncate);
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

		-- output: selection sums, shift and round
		PROC_STAGE_D: process(t3_r, i_s0, i_s1)
			variable v0, v1, v2, v3, v4 : t_acc;
		begin
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

			o_X0.re <= f_exit_round(v0, i_s0);
			o_X1.re <= f_exit_round(v1, i_s0);
			o_X2.re <= f_exit_round(v2, i_s0);
			o_X3.re <= f_exit_round(v3, i_s0);
			o_X4.re <= f_exit_round(v4, i_s0);

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

			o_X0.im <= f_exit_round(v0, i_s0);
			o_X1.im <= f_exit_round(v1, i_s0);
			o_X2.im <= f_exit_round(v2, i_s0);
			o_X3.im <= f_exit_round(v3, i_s0);
			o_X4.im <= f_exit_round(v4, i_s0);
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
		-- cut 1: after the input adder level
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

		-- cut 2: t1 sums and the biased raw k6 product
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

		-- cut 3: trim, multiply by j and t2 adds
		PROC_STAGE_C: process(i_clk)
			variable t1_2 : t_cmplx_prod;
			variable t2 : t_cmplx_prod_array;
		begin
			if rising_edge(i_clk) then
				t1_2.im := resize(r1_r.re, t1_2.im, fixed_wrap, fixed_truncate);
				t1_2.re := resize(r1_r.im, t1_2.re, fixed_wrap, fixed_truncate);
				t1_2.re := resize(-t1_2.re, t1_2.re, fixed_wrap, fixed_truncate);

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

		-- output: selection, shift and round
		PROC_STAGE_D: process(t2_r, i_s0)
			variable v0, v1, v2 : t_acc;
		begin
			v0 := resize(t2_r(0).re, v0);
			case i_s0(0) is
				when '0'    => v1 := resize(t2_r(3).re, v1);
				when others => v1 := resize(t2_r(1).re, v1);
			end case;
			v2 := resize(t2_r(2).re, v2);
			o_X0.re <= f_exit_round(v0, '0' & i_s0(0));
			o_X1.re <= f_exit_round(v1, '0' & i_s0(0));
			o_X2.re <= f_exit_round(v2, '0' & i_s0(0));

			v0 := resize(t2_r(0).im, v0);
			case i_s0(0) is
				when '0'    => v1 := resize(t2_r(3).im, v1);
				when others => v1 := resize(t2_r(1).im, v1);
			end case;
			v2 := resize(t2_r(2).im, v2);
			o_X0.im <= f_exit_round(v0, '0' & i_s0(0));
			o_X1.im <= f_exit_round(v1, '0' & i_s0(0));
			o_X2.im <= f_exit_round(v2, '0' & i_s0(0));
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
				o_X0.re <= f_exit_round(v0, "00");
				o_X1.re <= f_exit_round(v1, "00");

				v0 := resize(t1(0).im, v0);
				v1 := resize(t1(1).im, v1);
				o_X0.im <= f_exit_round(v0, "00");
				o_X1.im <= f_exit_round(v1, "00");
			end if;
		end process PROC_STAGE_A;
	end generate GEN_2_PIPE;

end architecture rtl;
