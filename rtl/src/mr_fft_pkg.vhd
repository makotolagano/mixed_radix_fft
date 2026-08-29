library ieee;
use ieee.std_logic_1164.all;
use ieee.fixed_pkg.all;
use ieee.fixed_float_types.all;
use ieee.math_real.all;


-- Mixed Rardix FFT Package Declaration Section
package mr_fft_pkg is

	-- "18 bits in memory, wide in flight" (docs/datapath_width_convention.md):
	-- every stored/streamed word (FIFOs, BRAM, inter-stage links, DSP inputs)
	-- is the s2.16 data word; inside the preadder the datapath is wider and
	-- NEVER re-quantized mid-flight -- the only roundings are the half-up
	-- product trim (prod word) and the single half-up shift+round at the
	-- preadder exit (plus the rotator's one exit rounding). Matches the
	-- Python model's exit_round mode (MixedRadix_PreAdder_FXP).
	constant c_fxp_int_width 		 	 : integer := 2;
	constant c_fxp_frac_width 		 : integer := 16;
	constant c_fxp_word_width 		 : integer := c_fxp_int_width + c_fxp_frac_width;
	constant c_guard_bits 				 : integer := 3;
	-- wide word (preadder adder tree, pre-multiplier): exact for ANY s2.16
	-- inputs -- tree magnitudes stay below 10 (int 5 = +-16) and the s0
	-- algorithmic shift (>>1/>>2) needs 2 extra fraction bits (18)
	constant c_fxp_int_wide_width  : integer := 5;
	constant c_fxp_frac_wide_width : integer := 18;
	-- product word (trimmed coefficient products, t3 adds, output sums):
	-- products trimmed half-up to 22 fraction bits; sums reach ~22.3 (int 6)
	constant c_fxp_prod_int_width  : integer := 6;
	constant c_fxp_prod_frac_width : integer := 22;
	constant c_coeff_int_width 		 : integer := 2;
	constant c_coeff_frac_width 	 : integer := 16;
	constant c_twiddle_int_width   : integer := 2;
	constant c_twiddle_frac_width  : integer := 16;
	-- Datapath quantization convention, used by every resize/to_sfixed in the
	-- arithmetic operators, the preadder and the rotator. Must match the
	-- Python model's rounding (round-to-nearest, see model fxp conventions);
	-- change both together or bit-exactness breaks. (The twiddle ROM's octant
	-- negations keep their own saturate/truncate -- exact-code negation, not
	-- data quantization.)
	constant c_fxp_round : fixed_round_style_type := fixed_round;

	-- a stage's supported FFT configs: (radix, stage_FFT_size_N) per config
	type t_config is record
			radix : natural; -- 2, 3 or 5
			size  : natural; -- twiddle size N at this stage for this config
	end record;
	type t_config_arr is array (natural range <>) of t_config;

	type t_nat_arr is array (natural range <>) of natural;

	-- Supported FFT lengths are N = c_fft_size_base * 2**i * 3**j * 5**k
	-- <= c_max_fft_size (same predicate as model/decomposition_configs.py).
	-- Config count, slot count, per-slot capability and per-slot config
	-- tables (see mr_fft_cfg_pkg) are all derived from these two numbers at
	-- elaboration time -- nothing else is hand-written.
	constant c_fft_size_base : natural := 12;
	constant c_max_fft_size  : natural := 3300;

	-- number of supported FFT lengths
	function f_num_configs return natural;
	-- supported FFT lengths, ascending; index = config index (i_config_sel)
	function f_fft_sizes return t_nat_arr;
	-- final-scaler constants, index-aligned with f_fft_sizes: code =
	-- round(2**c_scale_frac * 2**total_shift / N), the residual gain that
	-- turns the chain output X*2**(-total_shift) into the X/N convention.
	-- Under the FIXED schedule c = (4/3)**#3s * (8/5)**#5s, ranging from 4/3
	-- up to ~5.7 -- hence the s4.21 format (4 int bits incl sign, +-8). At
	-- 21 fraction bits the constant's relative error is ~-135 dB, far below
	-- every measured config (an 18-bit constant would sit ABOVE the best
	-- ones). The 25-bit code rides the DSP A port; the data word takes the
	-- 18-bit B port. Round-to-nearest is unambiguous: every supported N has
	-- a factor 3, so 2**x/N is never an exact half-integer.
	constant c_scale_int  : natural := 4;
	constant c_scale_frac : natural := 21;
	function f_fft_scales return t_nat_arr;
	-- minimal slot counts for the mid-bypass pipeline layout
	-- (docs/reconfigurable_fft_midbypass_architecture.md section 2):
	-- radix235 = max(#5s), radix23 = max(#3s+#5s) - radix235,
	-- radix2 = max(#2s+#3s+#5s) - radix23 - radix235
	function f_num_r2_slots   return natural;
	function f_num_r23_slots  return natural;
	function f_num_r235_slots return natural;
	function f_num_slots      return natural;
	-- capability of a pipeline slot (G_CAPABILITY encoding: 0=radix2,
	-- 1=radix23, 2=radix235); radix2 slots first, then radix23, then radix235
	function f_slot_capability(slot : natural) return natural;
	function f_ram_style(depth : natural) return string;
	-- (radix, size) of every config at one pipeline slot; bypassed configs
	-- get (radix=>1, size=>1). Use as G_CONFIGS of the stage at that slot.
	function f_slot_configs(slot : natural) return t_config_arr;

	type t_cmplx is record
		re : sfixed(c_fxp_int_width-1 downto -c_fxp_frac_width);
		im : sfixed(c_fxp_int_width-1 downto -c_fxp_frac_width);
  end record t_cmplx;
 
	type t_cmplx_wide is record
		re : sfixed(c_fxp_int_wide_width-1 downto -c_fxp_frac_wide_width);
		im : sfixed(c_fxp_int_wide_width-1 downto -c_fxp_frac_wide_width);
  end record t_cmplx_wide;

	type t_cmplx_coeff is record
		re : sfixed(c_coeff_int_width-1 downto -c_coeff_frac_width);
		im : sfixed(c_coeff_int_width-1 downto -c_coeff_frac_width);
	end record t_cmplx_coeff;

	type t_cmplx_prod is record
		re : sfixed(c_fxp_prod_int_width-1 downto -c_fxp_prod_frac_width);
		im : sfixed(c_fxp_prod_int_width-1 downto -c_fxp_prod_frac_width);
	end record t_cmplx_prod;

	type t_cmplx_mult is record
		re : sfixed(c_fxp_int_width+c_twiddle_int_width+1-1 downto -(c_fxp_frac_width+c_twiddle_frac_width));
		im : sfixed(c_fxp_int_width+c_twiddle_int_width+1-1 downto -(c_fxp_frac_width+c_twiddle_frac_width));
  end record t_cmplx_mult;

  type t_cmplx_wide_mult is record
		re : sfixed(c_fxp_int_wide_width+c_coeff_int_width+1-1 downto -(c_fxp_frac_wide_width+c_coeff_frac_width));
		im : sfixed(c_fxp_int_wide_width+c_coeff_int_width+1-1 downto -(c_fxp_frac_wide_width+c_coeff_frac_width));
  end record t_cmplx_wide_mult;

  type t_cmplx_twiddle is record
		re : sfixed(c_twiddle_int_width-1 downto -c_twiddle_frac_width);
		im : sfixed(c_twiddle_int_width-1 downto -c_twiddle_frac_width);
  end record t_cmplx_twiddle;

	-- Declare arithmetic operator prototypes for t_cmplx so they are
	-- visible at analysis time to units that `use` this package.
	function clogb2(n : integer) return integer;
	function get_max_radix(capability : natural) return natural;
	function get_delay_cnt(configs : t_config_arr) return natural;
	function get_max_size(configs : t_config_arr) return natural;
	-- preadder pipeline depth per capability (0 when not pipelined); the
	-- stage's valid pipeline and the preadder register cuts both use this
	function preadder_latency(capability : natural; pipelined : boolean) return natural;
	-- rotator pipeline depth (0 when combinational); the stage's skid sizing
	-- and the rotator register cuts both use this
	function rotator_latency(pipelined : boolean) return natural;
	function "+" (left, right : t_cmplx_wide) return t_cmplx_wide;
	function "-" (left, right : t_cmplx_wide) return t_cmplx_wide;
	function "*" (left : t_cmplx_wide; right : t_cmplx_coeff) return t_cmplx_wide;
	function "*" (left : t_cmplx; right : t_cmplx_twiddle) return t_cmplx;
	function resize(arg : t_cmplx; size_res : t_cmplx_wide) return t_cmplx_wide;
	function resize(arg : t_cmplx_wide; size_res : t_cmplx) return t_cmplx;
	function shift_right(arg : t_cmplx_wide; shift_amount : integer) return t_cmplx_wide;

	constant c_k2_re : sfixed(c_coeff_int_width-1 downto -c_coeff_frac_width):=
		to_sfixed(0.5 * (COS(MATH_2_PI / 5.0) - COS(2.0 * MATH_2_PI / 5.0)), c_coeff_int_width-1, -c_coeff_frac_width, fixed_wrap, c_fxp_round);

	constant c_k2 : t_cmplx_coeff := (
		re => c_k2_re,
		im => (others => '0')
	);

	constant c_k3_im : sfixed(c_coeff_int_width-1 downto -c_coeff_frac_width) :=
		to_sfixed(SIN(2.0 * MATH_2_PI / 5.0) - SIN(MATH_2_PI / 5.0), c_coeff_int_width-1, -c_coeff_frac_width, fixed_wrap, c_fxp_round);

	constant c_k3 : t_cmplx_coeff := (
		re => (others => '0'),
		im => c_k3_im
	);

	constant c_k4_im : sfixed(c_coeff_int_width-1 downto -c_coeff_frac_width) :=
		to_sfixed(-SIN(2.0 * MATH_2_PI / 5.0), c_coeff_int_width-1, -c_coeff_frac_width, fixed_wrap, c_fxp_round);

	constant c_k4 : t_cmplx_coeff := (
		re => (others => '0'),
		im => c_k4_im
	);

	constant c_k5_im : sfixed(c_coeff_int_width-1 downto -c_coeff_frac_width) :=
		to_sfixed(SIN(2.0 * MATH_2_PI / 5.0) + SIN(MATH_2_PI / 5.0), c_coeff_int_width-1, -c_coeff_frac_width, fixed_wrap, c_fxp_round);

	constant c_k5 : t_cmplx_coeff := (
		re => (others => '0'),
		im => c_k5_im
	);

	constant c_k6_re : sfixed(c_coeff_int_width-1 downto -c_coeff_frac_width) :=
		to_sfixed(-SQRT(3.0) / 2.0, c_coeff_int_width-1, -c_coeff_frac_width, fixed_wrap, c_fxp_round);

	constant c_k6 : t_cmplx_coeff := (
		re => c_k6_re,
		im => (others => '0')
	);


end package mr_fft_pkg;
 
-- Mixed Radix FFT Package Body Section
package body mr_fft_pkg is

	-- ceil(log2(n)), minimum of 1 bit.
	function clogb2(n : integer) return integer is
		variable res : integer := 0;
		variable v   : integer := n - 1;
	begin
		while v > 0 loop
			res := res + 1;
			v   := v / 2;
		end loop;
		if res = 0 then
			res := 1;
		end if;
		return res;
	end function;

  function get_max_radix(capability : natural) return natural is
  begin
    if capability = 2 then
      return 5;
		elsif capability = 1 then
      return 3;
    elsif capability = 0 then
      return 2;
    else
      return 0;
    end if;
  end function get_max_radix;

	function get_delay_cnt(configs : t_config_arr) return natural is
		variable max_N : natural := 0;
		variable radix : natural := 0;
	begin
		for i in configs'range loop
			if configs(i).size > max_N then
				max_N := configs(i).size;
				radix := configs(i).radix;
			end if;
		end loop;
		return (max_N/radix);
	end function get_delay_cnt;

	function get_max_size(configs : t_config_arr) return natural is
		variable m : natural := 1;
	begin
		for i in configs'range loop
			if configs(i).size > m then
				m := configs(i).size;
			end if;
		end loop;
		return m;
	end function get_max_size;

	function preadder_latency(capability : natural; pipelined : boolean) return natural is
	begin
		if not pipelined then
			return 0;
		end if;
		case capability is
			when 2      => return 3;  -- radix235: t1 adds / raw products / round + t3 adds
			when 1      => return 3;  -- radix23:  t0 adds / raw k6 product / round + t2 adds
			when others => return 1;  -- radix2:   registered butterfly outputs
		end case;
	end function preadder_latency;

	function rotator_latency(pipelined : boolean) return natural is
	begin
		if pipelined then
			return 5;  -- operand regs x2, M regs, post-adder/P regs, rounding reg
		else
			return 0;
		end if;
	end function rotator_latency;

	-- number of factors p in n
	function f_count_factor(n : natural; p : natural) return natural is
		variable v : natural := n;
		variable c : natural := 0;
	begin
		while v mod p = 0 loop
			v := v / p;
			c := c + 1;
		end loop;
		return c;
	end function f_count_factor;

	function f_num_configs return natural is
		variable n2, n3, n5 : natural;
		variable cnt        : natural := 0;
	begin
		n2 := c_fft_size_base;
		while n2 <= c_max_fft_size loop
			n3 := n2;
			while n3 <= c_max_fft_size loop
				n5 := n3;
				while n5 <= c_max_fft_size loop
					cnt := cnt + 1;
					n5  := n5 * 5;
				end loop;
				n3 := n3 * 3;
			end loop;
			n2 := n2 * 2;
		end loop;
		return cnt;
	end function f_num_configs;

	function f_fft_sizes return t_nat_arr is
		variable sizes      : t_nat_arr(0 to f_num_configs - 1);
		variable idx        : natural := 0;
		variable n2, n3, n5 : natural;
		variable tmp        : natural;
		variable j          : natural;
	begin
		n2 := c_fft_size_base;
		while n2 <= c_max_fft_size loop
			n3 := n2;
			while n3 <= c_max_fft_size loop
				n5 := n3;
				while n5 <= c_max_fft_size loop
					sizes(idx) := n5;
					idx        := idx + 1;
					n5         := n5 * 5;
				end loop;
				n3 := n3 * 3;
			end loop;
			n2 := n2 * 2;
		end loop;
		-- insertion sort ascending: config index = rank of N (model ordering)
		for i in 1 to sizes'high loop
			tmp := sizes(i);
			j   := i;
			while j > 0 loop
				exit when sizes(j - 1) <= tmp;
				sizes(j) := sizes(j - 1);
				j        := j - 1;
			end loop;
			sizes(j) := tmp;
		end loop;
		return sizes;
	end function f_fft_sizes;

	function f_fft_scales return t_nat_arr is
		constant sizes : t_nat_arr := f_fft_sizes;
		variable t     : t_nat_arr(sizes'range);
		variable sh    : natural;
	begin
		for i in sizes'range loop
			-- fixed scaling schedule: 1 bit per radix-2, 2 per radix-3, 3 per radix-5
			sh := f_count_factor(sizes(i), 2) + 2 * f_count_factor(sizes(i), 3)
			      + 3 * f_count_factor(sizes(i), 5);
			t(i) := natural(round(2.0 ** c_scale_frac * 2.0 ** sh / real(sizes(i))));
		end loop;
		return t;
	end function f_fft_scales;

	function f_num_r235_slots return natural is
		constant sizes : t_nat_arr := f_fft_sizes;
		variable m     : natural   := 0;
	begin
		for i in sizes'range loop
			if f_count_factor(sizes(i), 5) > m then
				m := f_count_factor(sizes(i), 5);
			end if;
		end loop;
		return m;
	end function f_num_r235_slots;

	function f_num_r23_slots return natural is
		constant sizes : t_nat_arr := f_fft_sizes;
		variable f     : natural;
		variable m     : natural := 0;
	begin
		for i in sizes'range loop
			f := f_count_factor(sizes(i), 3) + f_count_factor(sizes(i), 5);
			if f > m then
				m := f;
			end if;
		end loop;
		return m - f_num_r235_slots;
	end function f_num_r23_slots;

	function f_num_r2_slots return natural is
		constant sizes : t_nat_arr := f_fft_sizes;
		variable f     : natural;
		variable m     : natural := 0;
	begin
		for i in sizes'range loop
			f := f_count_factor(sizes(i), 2) + f_count_factor(sizes(i), 3) + f_count_factor(sizes(i), 5);
			if f > m then
				m := f;
			end if;
		end loop;
		return m - f_num_r23_slots - f_num_r235_slots;
	end function f_num_r2_slots;

	function f_num_slots return natural is
	begin
		return f_num_r2_slots + f_num_r23_slots + f_num_r235_slots;
	end function f_num_slots;

	function f_slot_capability(slot : natural) return natural is
	begin
		if slot < f_num_r2_slots then
			return 0;
		elsif slot < f_num_r2_slots + f_num_r23_slots then
			return 1;
		else
			return 2;
		end if;
	end function f_slot_capability;

	function f_ram_style(depth : natural) return string is
	begin
		if depth >= 128 then
			return "block";
		else
			return "distributed";
		end if;
	end function f_ram_style;

	-- Mid-bypass placement: each config's ascending radix list [2..][3..][5..]
	-- goes to the earliest slot whose capability admits the radix (2: any
	-- slot, 3: radix23/radix235, 5: radix235 only); remaining slots are
	-- bypassed. size at a used slot is the remaining DIF sub-FFT length: the
	-- product of the radices placed at this and all later slots.
	function f_slot_configs(slot : natural) return t_config_arr is
		constant sizes      : t_nat_arr := f_fft_sizes;
		constant num_slots  : natural   := f_num_slots;
		constant first_r23  : natural   := f_num_r2_slots;
		constant first_r235 : natural   := f_num_r2_slots + f_num_r23_slots;
		variable radix_at   : t_nat_arr(0 to num_slots - 1);
		variable pos        : natural;
		variable prod       : natural;
		variable cfgs       : t_config_arr(0 to sizes'length - 1);
	begin
		assert slot < num_slots
			report "f_slot_configs: slot " & integer'image(slot) & " out of range"
			severity failure;
		for c in sizes'range loop
			radix_at := (others => 1);
			pos      := 0;
			for r in 1 to f_count_factor(sizes(c), 2) loop
				radix_at(pos) := 2;
				pos           := pos + 1;
			end loop;
			if pos < first_r23 then
				pos := first_r23;
			end if;
			for r in 1 to f_count_factor(sizes(c), 3) loop
				radix_at(pos) := 3;
				pos           := pos + 1;
			end loop;
			if pos < first_r235 then
				pos := first_r235;
			end if;
			for r in 1 to f_count_factor(sizes(c), 5) loop
				radix_at(pos) := 5;
				pos           := pos + 1;
			end loop;
			-- remaining sub-FFT length seen at `slot` (bypass slots contribute 1)
			prod := 1;
			for s in num_slots - 1 downto slot loop
				prod := prod * radix_at(s);
			end loop;
			if radix_at(slot) = 1 then
				cfgs(c) := (radix => 1, size => 1);
			else
				cfgs(c) := (radix => radix_at(slot), size => prod);
			end if;
			-- self-check: the placed radices must multiply back to N
			for s in slot - 1 downto 0 loop
				prod := prod * radix_at(s);
			end loop;
			assert prod = sizes(c)
				report "f_slot_configs: placement does not decompose N=" & integer'image(sizes(c))
				severity failure;
		end loop;
		return cfgs;
	end function f_slot_configs;

	function "+" (left, right : t_cmplx_wide) return t_cmplx_wide is
    variable result : t_cmplx_wide;
	begin
		result.re := resize(left.re + right.re, c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, c_fxp_round);
		result.im := resize(left.im + right.im, c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, c_fxp_round);
		return result;
	end function;

	function "-" (left, right : t_cmplx_wide) return t_cmplx_wide is
		variable result : t_cmplx_wide;
	begin
		result.re := resize(left.re - right.re, c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, c_fxp_round);
		result.im := resize(left.im - right.im, c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, c_fxp_round);
		return result;
	end function;

	function "*" (left : t_cmplx_wide; right : t_cmplx_coeff) return t_cmplx_wide is
		variable result 		 : t_cmplx_wide;
		variable mult_result : t_cmplx_wide_mult;
	begin
		mult_result.re := left.re * right.re - left.im * right.im;
		mult_result.im := left.re * right.im + left.im * right.re;
		result.re := resize(mult_result.re, result.re, fixed_wrap, c_fxp_round);
		result.im := resize(mult_result.im, result.im, fixed_wrap, c_fxp_round);
		return result;
	end function;

	-- rotator product: full-precision 4-multiplier complex product with ONE
	-- half-up exit rounding (+half LSB then truncate, wrap) back to the
	-- memory word -- matches the model's rounding='half_up' rotator quantizer
	function "*" (left : t_cmplx; right : t_cmplx_twiddle) return t_cmplx is
		variable result 		 : t_cmplx;
		variable mult_result : t_cmplx_mult;
		constant c_half : sfixed(0 downto -(c_fxp_frac_width + 1)) :=
			to_sfixed(2.0 ** (-(c_fxp_frac_width + 1)), 0, -(c_fxp_frac_width + 1));
	begin
		mult_result.re := left.re * right.re - left.im * right.im;
		mult_result.im := left.re * right.im + left.im * right.re;
		result.re := resize(mult_result.re + c_half, result.re, fixed_wrap, fixed_truncate);
		result.im := resize(mult_result.im + c_half, result.im, fixed_wrap, fixed_truncate);
		return result;
	end function;

	function resize(arg : t_cmplx; size_res : t_cmplx_wide) return t_cmplx_wide is
		variable result : t_cmplx_wide;
	begin
		result.re := resize(arg.re, result.re, fixed_wrap, c_fxp_round);
		result.im := resize(arg.im, result.im, fixed_wrap, c_fxp_round);
		return result;
	end function;

	function resize(arg : t_cmplx_wide; size_res : t_cmplx) return t_cmplx is
		variable result : t_cmplx;
	begin
		result.re := resize(arg.re, result.re, fixed_wrap, c_fxp_round);
		result.im := resize(arg.im, result.im, fixed_wrap, c_fxp_round);
		return result;
	end function;

	function shift_right(arg : t_cmplx_wide; shift_amount : integer) return t_cmplx_wide is
		variable result : t_cmplx_wide;
	begin
		result.re := resize(shift_right(arg.re, shift_amount), result.re, fixed_wrap, c_fxp_round);
		result.im := resize(shift_right(arg.im, shift_amount), result.im, fixed_wrap, c_fxp_round);
		return result;
	end function;
 
end package body mr_fft_pkg;