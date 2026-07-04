library ieee;
use ieee.std_logic_1164.all;
use ieee.fixed_pkg.all;
use ieee.fixed_float_types.all;
use ieee.math_real.all;


-- Mixed Rardix FFT Package Declaration Section
package mr_fft_pkg is

	constant c_data_width : integer := 16;

	constant c_fxp_int_width 		 	 : integer := 5;
	constant c_fxp_frac_width 		 : integer := 12;
	constant c_guard_bits 				 : integer := 3;
	constant c_fxp_int_wide_width  : integer := c_fxp_int_width + c_guard_bits;
	constant c_fxp_frac_wide_width : integer := c_fxp_frac_width + c_guard_bits;
	
	type t_cmplx is record
		re : sfixed(c_fxp_int_width-1 downto -c_fxp_frac_width);
		im : sfixed(c_fxp_int_width-1 downto -c_fxp_frac_width);
  end record t_cmplx;
 
	type t_cmplx_wide is record
		re : sfixed(c_fxp_int_wide_width-1 downto -c_fxp_frac_wide_width);
		im : sfixed(c_fxp_int_wide_width-1 downto -c_fxp_frac_wide_width);
  end record t_cmplx_wide;
  
	type t_cmplx_wide_mult is record
		re : sfixed(c_fxp_int_wide_width*2+1-1 downto -c_fxp_frac_wide_width*2);
		im : sfixed(c_fxp_int_wide_width*2+1-1 downto -c_fxp_frac_wide_width*2);
  end record t_cmplx_wide_mult;

	-- Declare arithmetic operator prototypes for t_cmplx so they are
	-- visible at analysis time to units that `use` this package.
	function "+" (left, right : t_cmplx_wide) return t_cmplx_wide;
	function "-" (left, right : t_cmplx_wide) return t_cmplx_wide;
	function "*" (left, right : t_cmplx_wide) return t_cmplx_wide;
	function resize(arg : t_cmplx; size_res : t_cmplx_wide) return t_cmplx_wide;
	function resize(arg : t_cmplx_wide; size_res : t_cmplx) return t_cmplx;
	function shift_right(arg : t_cmplx_wide; shift_amount : integer) return t_cmplx_wide;

	constant c_k2_re : sfixed(c_fxp_int_wide_width-1 downto -c_fxp_frac_wide_width):=
		to_sfixed(0.5 * (COS(MATH_2_PI / 5.0) - COS(2.0 * MATH_2_PI / 5.0)), c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, fixed_truncate);

	constant c_k2 : t_cmplx_wide := (
		re => c_k2_re,
		im => (others => '0')
	);

	constant c_k3_im : sfixed(c_fxp_int_wide_width-1 downto -c_fxp_frac_wide_width) :=
		to_sfixed(SIN(2.0 * MATH_2_PI / 5.0) - SIN(MATH_2_PI / 5.0), c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, fixed_truncate);

	constant c_k3 : t_cmplx_wide := (
		re => (others => '0'),
		im => c_k3_im
	);

	constant c_k4_im : sfixed(c_fxp_int_wide_width-1 downto -c_fxp_frac_wide_width) :=
		to_sfixed(-SIN(2.0 * MATH_2_PI / 5.0), c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, fixed_truncate);

	constant c_k4 : t_cmplx_wide := (
		re => (others => '0'),
		im => c_k4_im
	);

	constant c_k5_im : sfixed(c_fxp_int_wide_width-1 downto -c_fxp_frac_wide_width) :=
		to_sfixed(SIN(2.0 * MATH_2_PI / 5.0) + SIN(MATH_2_PI / 5.0), c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, fixed_truncate);

	constant c_k5 : t_cmplx_wide := (
		re => (others => '0'),
		im => c_k5_im
	);

	constant c_k6_re : sfixed(c_fxp_int_wide_width-1 downto -c_fxp_frac_wide_width) :=
		to_sfixed(-SQRT(3.0) / 2.0, c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, fixed_truncate);

	constant c_k6 : t_cmplx_wide := (
		re => c_k6_re,
		im => (others => '0')
	);

end package mr_fft_pkg;
 
-- Mixed Radix FFT Package Body Section
package body mr_fft_pkg is

	function "+" (left, right : t_cmplx_wide) return t_cmplx_wide is
    variable result : t_cmplx_wide;
	begin
		result.re := resize(left.re + right.re, c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, fixed_truncate);
		result.im := resize(left.im + right.im, c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, fixed_truncate);
		return result;
	end function;

	function "-" (left, right : t_cmplx_wide) return t_cmplx_wide is
		variable result : t_cmplx_wide;
	begin
		result.re := resize(left.re - right.re, c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, fixed_truncate);
		result.im := resize(left.im - right.im, c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, fixed_truncate);
		return result;
	end function;

	function "*" (left, right : t_cmplx_wide) return t_cmplx_wide is
		variable result 		 : t_cmplx_wide;
		variable mult_result : t_cmplx_wide_mult;
	begin
		mult_result.re := left.re * right.re - left.im * right.im;
		mult_result.im := left.re * right.im + left.im * right.re;
		result.re := resize(mult_result.re, result.re, fixed_wrap, fixed_truncate);
		result.im := resize(mult_result.im, result.im, fixed_wrap, fixed_truncate);
		return result;
	end function;

	function resize(arg : t_cmplx; size_res : t_cmplx_wide) return t_cmplx_wide is
		variable result : t_cmplx_wide;
	begin
		result.re := resize(arg.re, result.re, fixed_wrap, fixed_truncate);
		result.im := resize(arg.im, result.im, fixed_wrap, fixed_truncate);
		return result;
	end function;

	function resize(arg : t_cmplx_wide; size_res : t_cmplx) return t_cmplx is
		variable result : t_cmplx;
	begin
		result.re := resize(arg.re, result.re, fixed_wrap, fixed_truncate);
		result.im := resize(arg.im, result.im, fixed_wrap, fixed_truncate);
		return result;
	end function;

	function shift_right(arg : t_cmplx_wide; shift_amount : integer) return t_cmplx_wide is
		variable result : t_cmplx_wide;
	begin
		result.re := resize(shift_right(arg.re, shift_amount), result.re, fixed_wrap, fixed_truncate);
		result.im := resize(shift_right(arg.im, shift_amount), result.im, fixed_wrap, fixed_truncate);
		return result;
	end function;
 
end package body mr_fft_pkg;