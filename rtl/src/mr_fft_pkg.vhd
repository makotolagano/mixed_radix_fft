library ieee;
use ieee.std_logic_1164.all;
use ieee.fixed_pkg.all;
use ieee.fixed_float_types.all;
use ieee.math_real.all;


-- Mixed Rardix FFT Package Declaration Section
package mr_fft_pkg is

	constant c_fxp_int_width 		 	 : integer := 5;
	constant c_fxp_frac_width 		 : integer := 13;
	constant c_fxp_word_width 		 : integer := c_fxp_int_width + c_fxp_frac_width;
	constant c_guard_bits 				 : integer := 3;
	-- The Python model (MixedRadix_PreAdder_FXP) quantizes ports and every
	-- intermediate at ONE dtype (inner_type); the wide format must therefore
	-- equal the data format to stay bit-exact with it.
	constant c_fxp_int_wide_width  : integer := c_fxp_int_width;
	constant c_fxp_frac_wide_width : integer := c_fxp_frac_width;
	constant c_coeff_int_width 		 : integer := 2;
	constant c_coeff_frac_width 	 : integer := 16;
	constant c_twiddle_int_width   : integer := 2;
	constant c_twiddle_frac_width  : integer := 16;

	-- a stage's supported FFT configs: (radix, stage_FFT_size_N) per config
	type t_config is record
			radix : natural; -- 2, 3 or 5
			size  : natural; -- twiddle size N at this stage for this config
	end record;
	type t_config_arr is array (natural range <>) of t_config;
	
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
	function "+" (left, right : t_cmplx_wide) return t_cmplx_wide;
	function "-" (left, right : t_cmplx_wide) return t_cmplx_wide;
	function "*" (left : t_cmplx_wide; right : t_cmplx_coeff) return t_cmplx_wide;
	function "*" (left : t_cmplx; right : t_cmplx_twiddle) return t_cmplx;
	function resize(arg : t_cmplx; size_res : t_cmplx_wide) return t_cmplx_wide;
	function resize(arg : t_cmplx_wide; size_res : t_cmplx) return t_cmplx;
	function shift_right(arg : t_cmplx_wide; shift_amount : integer) return t_cmplx_wide;

	constant c_k2_re : sfixed(c_coeff_int_width-1 downto -c_coeff_frac_width):=
		to_sfixed(0.5 * (COS(MATH_2_PI / 5.0) - COS(2.0 * MATH_2_PI / 5.0)), c_coeff_int_width-1, -c_coeff_frac_width, fixed_wrap, fixed_round);

	constant c_k2 : t_cmplx_coeff := (
		re => c_k2_re,
		im => (others => '0')
	);

	constant c_k3_im : sfixed(c_coeff_int_width-1 downto -c_coeff_frac_width) :=
		to_sfixed(SIN(2.0 * MATH_2_PI / 5.0) - SIN(MATH_2_PI / 5.0), c_coeff_int_width-1, -c_coeff_frac_width, fixed_wrap, fixed_round);

	constant c_k3 : t_cmplx_coeff := (
		re => (others => '0'),
		im => c_k3_im
	);

	constant c_k4_im : sfixed(c_coeff_int_width-1 downto -c_coeff_frac_width) :=
		to_sfixed(-SIN(2.0 * MATH_2_PI / 5.0), c_coeff_int_width-1, -c_coeff_frac_width, fixed_wrap, fixed_round);

	constant c_k4 : t_cmplx_coeff := (
		re => (others => '0'),
		im => c_k4_im
	);

	constant c_k5_im : sfixed(c_coeff_int_width-1 downto -c_coeff_frac_width) :=
		to_sfixed(SIN(2.0 * MATH_2_PI / 5.0) + SIN(MATH_2_PI / 5.0), c_coeff_int_width-1, -c_coeff_frac_width, fixed_wrap, fixed_round);

	constant c_k5 : t_cmplx_coeff := (
		re => (others => '0'),
		im => c_k5_im
	);

	constant c_k6_re : sfixed(c_coeff_int_width-1 downto -c_coeff_frac_width) :=
		to_sfixed(-SQRT(3.0) / 2.0, c_coeff_int_width-1, -c_coeff_frac_width, fixed_wrap, fixed_round);

	constant c_k6 : t_cmplx_coeff := (
		re => c_k6_re,
		im => (others => '0')
	);

	-- Slot 0: radix235  (declared radix 2/3/5)
	constant c_stage0_cfgs : t_config_arr := (
		 0 => (radix => 3, size =>   12),  -- N=12   depth=4
		 1 => (radix => 3, size =>   24),  -- N=24   depth=8
		 2 => (radix => 3, size =>   36),  -- N=36   depth=12
		 3 => (radix => 3, size =>   48),  -- N=48   depth=16
		 4 => (radix => 5, size =>   60),  -- N=60   depth=12
		 5 => (radix => 3, size =>   72),  -- N=72   depth=24
		 6 => (radix => 3, size =>   96),  -- N=96   depth=32
		 7 => (radix => 3, size =>  108),  -- N=108  depth=36
		 8 => (radix => 5, size =>  120),  -- N=120  depth=24
		 9 => (radix => 3, size =>  144),  -- N=144  depth=48
		10 => (radix => 5, size =>  180),  -- N=180  depth=36
		11 => (radix => 3, size =>  192),  -- N=192  depth=64
		12 => (radix => 3, size =>  216),  -- N=216  depth=72
		13 => (radix => 5, size =>  240),  -- N=240  depth=48
		14 => (radix => 3, size =>  288),  -- N=288  depth=96
		15 => (radix => 5, size =>  300),  -- N=300  depth=60
		16 => (radix => 3, size =>  324),  -- N=324  depth=108
		17 => (radix => 5, size =>  360),  -- N=360  depth=72
		18 => (radix => 3, size =>  384),  -- N=384  depth=128
		19 => (radix => 3, size =>  432),  -- N=432  depth=144
		20 => (radix => 5, size =>  480),  -- N=480  depth=96
		21 => (radix => 5, size =>  540),  -- N=540  depth=108
		22 => (radix => 3, size =>  576),  -- N=576  depth=192
		23 => (radix => 5, size =>  600),  -- N=600  depth=120
		24 => (radix => 3, size =>  648),  -- N=648  depth=216
		25 => (radix => 5, size =>  720),  -- N=720  depth=144
		26 => (radix => 3, size =>  768),  -- N=768  depth=256
		27 => (radix => 3, size =>  864),  -- N=864  depth=288
		28 => (radix => 5, size =>  900),  -- N=900  depth=180
		29 => (radix => 5, size =>  960),  -- N=960  depth=192
		30 => (radix => 3, size =>  972),  -- N=972  depth=324
		31 => (radix => 5, size => 1080),  -- N=1080 depth=216
		32 => (radix => 3, size => 1152),  -- N=1152 depth=384
		33 => (radix => 5, size => 1200),  -- N=1200 depth=240
		34 => (radix => 3, size => 1296),  -- N=1296 depth=432
		35 => (radix => 5, size => 1440),  -- N=1440 depth=288
		36 => (radix => 5, size => 1500),  -- N=1500 depth=300
		37 => (radix => 3, size => 1536),  -- N=1536 depth=512
		38 => (radix => 5, size => 1620),  -- N=1620 depth=324
		39 => (radix => 3, size => 1728),  -- N=1728 depth=576
		40 => (radix => 5, size => 1800),  -- N=1800 depth=360
		41 => (radix => 5, size => 1920),  -- N=1920 depth=384
		42 => (radix => 3, size => 1944),  -- N=1944 depth=648
		43 => (radix => 5, size => 2160),  -- N=2160 depth=432
		44 => (radix => 3, size => 2304),  -- N=2304 depth=768
		45 => (radix => 5, size => 2400),  -- N=2400 depth=480
		46 => (radix => 3, size => 2592),  -- N=2592 depth=864
		47 => (radix => 5, size => 2700),  -- N=2700 depth=540
		48 => (radix => 5, size => 2880),  -- N=2880 depth=576
		49 => (radix => 3, size => 2916),  -- N=2916 depth=972
		50 => (radix => 5, size => 3000),  -- N=3000 depth=600
		51 => (radix => 3, size => 3072),  -- N=3072 depth=1024
		52 => (radix => 5, size => 3240)   -- N=3240 depth=648
	);

	-- Slot 1: radix235  (declared radix 2/3/5)
	constant c_stage1_cfgs : t_config_arr := (
		 0 => (radix => 2, size =>    4),  -- N=12   depth=2
		 1 => (radix => 2, size =>    8),  -- N=24   depth=4
		 2 => (radix => 3, size =>   12),  -- N=36   depth=4
		 3 => (radix => 2, size =>   16),  -- N=48   depth=8
		 4 => (radix => 3, size =>   12),  -- N=60   depth=4
		 5 => (radix => 3, size =>   24),  -- N=72   depth=8
		 6 => (radix => 2, size =>   32),  -- N=96   depth=16
		 7 => (radix => 3, size =>   36),  -- N=108  depth=12
		 8 => (radix => 3, size =>   24),  -- N=120  depth=8
		 9 => (radix => 3, size =>   48),  -- N=144  depth=16
		10 => (radix => 3, size =>   36),  -- N=180  depth=12
		11 => (radix => 2, size =>   64),  -- N=192  depth=32
		12 => (radix => 3, size =>   72),  -- N=216  depth=24
		13 => (radix => 3, size =>   48),  -- N=240  depth=16
		14 => (radix => 3, size =>   96),  -- N=288  depth=32
		15 => (radix => 5, size =>   60),  -- N=300  depth=12
		16 => (radix => 3, size =>  108),  -- N=324  depth=36
		17 => (radix => 3, size =>   72),  -- N=360  depth=24
		18 => (radix => 2, size =>  128),  -- N=384  depth=64
		19 => (radix => 3, size =>  144),  -- N=432  depth=48
		20 => (radix => 3, size =>   96),  -- N=480  depth=32
		21 => (radix => 3, size =>  108),  -- N=540  depth=36
		22 => (radix => 3, size =>  192),  -- N=576  depth=64
		23 => (radix => 5, size =>  120),  -- N=600  depth=24
		24 => (radix => 3, size =>  216),  -- N=648  depth=72
		25 => (radix => 3, size =>  144),  -- N=720  depth=48
		26 => (radix => 2, size =>  256),  -- N=768  depth=128
		27 => (radix => 3, size =>  288),  -- N=864  depth=96
		28 => (radix => 5, size =>  180),  -- N=900  depth=36
		29 => (radix => 3, size =>  192),  -- N=960  depth=64
		30 => (radix => 3, size =>  324),  -- N=972  depth=108
		31 => (radix => 3, size =>  216),  -- N=1080 depth=72
		32 => (radix => 3, size =>  384),  -- N=1152 depth=128
		33 => (radix => 5, size =>  240),  -- N=1200 depth=48
		34 => (radix => 3, size =>  432),  -- N=1296 depth=144
		35 => (radix => 3, size =>  288),  -- N=1440 depth=96
		36 => (radix => 5, size =>  300),  -- N=1500 depth=60
		37 => (radix => 2, size =>  512),  -- N=1536 depth=256
		38 => (radix => 3, size =>  324),  -- N=1620 depth=108
		39 => (radix => 3, size =>  576),  -- N=1728 depth=192
		40 => (radix => 5, size =>  360),  -- N=1800 depth=72
		41 => (radix => 3, size =>  384),  -- N=1920 depth=128
		42 => (radix => 3, size =>  648),  -- N=1944 depth=216
		43 => (radix => 3, size =>  432),  -- N=2160 depth=144
		44 => (radix => 3, size =>  768),  -- N=2304 depth=256
		45 => (radix => 5, size =>  480),  -- N=2400 depth=96
		46 => (radix => 3, size =>  864),  -- N=2592 depth=288
		47 => (radix => 5, size =>  540),  -- N=2700 depth=108
		48 => (radix => 3, size =>  576),  -- N=2880 depth=192
		49 => (radix => 3, size =>  972),  -- N=2916 depth=324
		50 => (radix => 5, size =>  600),  -- N=3000 depth=120
		51 => (radix => 2, size => 1024),  -- N=3072 depth=512
		52 => (radix => 3, size =>  648)   -- N=3240 depth=216
	);

	-- Slot 2: radix235  (declared radix 2/3/5)
	constant c_stage2_cfgs : t_config_arr := (
		 0 => (radix => 2, size =>    2),  -- N=12   depth=1
		 1 => (radix => 2, size =>    4),  -- N=24   depth=2
		 2 => (radix => 2, size =>    4),  -- N=36   depth=2
		 3 => (radix => 2, size =>    8),  -- N=48   depth=4
		 4 => (radix => 2, size =>    4),  -- N=60   depth=2
		 5 => (radix => 2, size =>    8),  -- N=72   depth=4
		 6 => (radix => 2, size =>   16),  -- N=96   depth=8
		 7 => (radix => 3, size =>   12),  -- N=108  depth=4
		 8 => (radix => 2, size =>    8),  -- N=120  depth=4
		 9 => (radix => 2, size =>   16),  -- N=144  depth=8
		10 => (radix => 3, size =>   12),  -- N=180  depth=4
		11 => (radix => 2, size =>   32),  -- N=192  depth=16
		12 => (radix => 3, size =>   24),  -- N=216  depth=8
		13 => (radix => 2, size =>   16),  -- N=240  depth=8
		14 => (radix => 2, size =>   32),  -- N=288  depth=16
		15 => (radix => 3, size =>   12),  -- N=300  depth=4
		16 => (radix => 3, size =>   36),  -- N=324  depth=12
		17 => (radix => 3, size =>   24),  -- N=360  depth=8
		18 => (radix => 2, size =>   64),  -- N=384  depth=32
		19 => (radix => 3, size =>   48),  -- N=432  depth=16
		20 => (radix => 2, size =>   32),  -- N=480  depth=16
		21 => (radix => 3, size =>   36),  -- N=540  depth=12
		22 => (radix => 2, size =>   64),  -- N=576  depth=32
		23 => (radix => 3, size =>   24),  -- N=600  depth=8
		24 => (radix => 3, size =>   72),  -- N=648  depth=24
		25 => (radix => 3, size =>   48),  -- N=720  depth=16
		26 => (radix => 2, size =>  128),  -- N=768  depth=64
		27 => (radix => 3, size =>   96),  -- N=864  depth=32
		28 => (radix => 3, size =>   36),  -- N=900  depth=12
		29 => (radix => 2, size =>   64),  -- N=960  depth=32
		30 => (radix => 3, size =>  108),  -- N=972  depth=36
		31 => (radix => 3, size =>   72),  -- N=1080 depth=24
		32 => (radix => 2, size =>  128),  -- N=1152 depth=64
		33 => (radix => 3, size =>   48),  -- N=1200 depth=16
		34 => (radix => 3, size =>  144),  -- N=1296 depth=48
		35 => (radix => 3, size =>   96),  -- N=1440 depth=32
		36 => (radix => 5, size =>   60),  -- N=1500 depth=12
		37 => (radix => 2, size =>  256),  -- N=1536 depth=128
		38 => (radix => 3, size =>  108),  -- N=1620 depth=36
		39 => (radix => 3, size =>  192),  -- N=1728 depth=64
		40 => (radix => 3, size =>   72),  -- N=1800 depth=24
		41 => (radix => 2, size =>  128),  -- N=1920 depth=64
		42 => (radix => 3, size =>  216),  -- N=1944 depth=72
		43 => (radix => 3, size =>  144),  -- N=2160 depth=48
		44 => (radix => 2, size =>  256),  -- N=2304 depth=128
		45 => (radix => 3, size =>   96),  -- N=2400 depth=32
		46 => (radix => 3, size =>  288),  -- N=2592 depth=96
		47 => (radix => 3, size =>  108),  -- N=2700 depth=36
		48 => (radix => 3, size =>  192),  -- N=2880 depth=64
		49 => (radix => 3, size =>  324),  -- N=2916 depth=108
		50 => (radix => 5, size =>  120),  -- N=3000 depth=24
		51 => (radix => 2, size =>  512),  -- N=3072 depth=256
		52 => (radix => 3, size =>  216)   -- N=3240 depth=72
	);

	-- Slot 3: radix23  (declared radix 2/3)
	constant c_stage3_cfgs : t_config_arr := (
		 0 => (radix => 1, size =>    1),  -- N=12   bypass
		 1 => (radix => 2, size =>    2),  -- N=24   depth=1
		 2 => (radix => 2, size =>    2),  -- N=36   depth=1
		 3 => (radix => 2, size =>    4),  -- N=48   depth=2
		 4 => (radix => 2, size =>    2),  -- N=60   depth=1
		 5 => (radix => 2, size =>    4),  -- N=72   depth=2
		 6 => (radix => 2, size =>    8),  -- N=96   depth=4
		 7 => (radix => 2, size =>    4),  -- N=108  depth=2
		 8 => (radix => 2, size =>    4),  -- N=120  depth=2
		 9 => (radix => 2, size =>    8),  -- N=144  depth=4
		10 => (radix => 2, size =>    4),  -- N=180  depth=2
		11 => (radix => 2, size =>   16),  -- N=192  depth=8
		12 => (radix => 2, size =>    8),  -- N=216  depth=4
		13 => (radix => 2, size =>    8),  -- N=240  depth=4
		14 => (radix => 2, size =>   16),  -- N=288  depth=8
		15 => (radix => 2, size =>    4),  -- N=300  depth=2
		16 => (radix => 3, size =>   12),  -- N=324  depth=4
		17 => (radix => 2, size =>    8),  -- N=360  depth=4
		18 => (radix => 2, size =>   32),  -- N=384  depth=16
		19 => (radix => 2, size =>   16),  -- N=432  depth=8
		20 => (radix => 2, size =>   16),  -- N=480  depth=8
		21 => (radix => 3, size =>   12),  -- N=540  depth=4
		22 => (radix => 2, size =>   32),  -- N=576  depth=16
		23 => (radix => 2, size =>    8),  -- N=600  depth=4
		24 => (radix => 3, size =>   24),  -- N=648  depth=8
		25 => (radix => 2, size =>   16),  -- N=720  depth=8
		26 => (radix => 2, size =>   64),  -- N=768  depth=32
		27 => (radix => 2, size =>   32),  -- N=864  depth=16
		28 => (radix => 3, size =>   12),  -- N=900  depth=4
		29 => (radix => 2, size =>   32),  -- N=960  depth=16
		30 => (radix => 3, size =>   36),  -- N=972  depth=12
		31 => (radix => 3, size =>   24),  -- N=1080 depth=8
		32 => (radix => 2, size =>   64),  -- N=1152 depth=32
		33 => (radix => 2, size =>   16),  -- N=1200 depth=8
		34 => (radix => 3, size =>   48),  -- N=1296 depth=16
		35 => (radix => 2, size =>   32),  -- N=1440 depth=16
		36 => (radix => 3, size =>   12),  -- N=1500 depth=4
		37 => (radix => 2, size =>  128),  -- N=1536 depth=64
		38 => (radix => 3, size =>   36),  -- N=1620 depth=12
		39 => (radix => 2, size =>   64),  -- N=1728 depth=32
		40 => (radix => 3, size =>   24),  -- N=1800 depth=8
		41 => (radix => 2, size =>   64),  -- N=1920 depth=32
		42 => (radix => 3, size =>   72),  -- N=1944 depth=24
		43 => (radix => 3, size =>   48),  -- N=2160 depth=16
		44 => (radix => 2, size =>  128),  -- N=2304 depth=64
		45 => (radix => 2, size =>   32),  -- N=2400 depth=16
		46 => (radix => 3, size =>   96),  -- N=2592 depth=32
		47 => (radix => 3, size =>   36),  -- N=2700 depth=12
		48 => (radix => 2, size =>   64),  -- N=2880 depth=32
		49 => (radix => 3, size =>  108),  -- N=2916 depth=36
		50 => (radix => 3, size =>   24),  -- N=3000 depth=8
		51 => (radix => 2, size =>  256),  -- N=3072 depth=128
		52 => (radix => 3, size =>   72)   -- N=3240 depth=24
	);

	-- Slot 4: radix23  (declared radix 2/3)
	constant c_stage4_cfgs : t_config_arr := (
		 0 => (radix => 1, size =>    1),  -- N=12   bypass
		 1 => (radix => 1, size =>    1),  -- N=24   bypass
		 2 => (radix => 1, size =>    1),  -- N=36   bypass
		 3 => (radix => 2, size =>    2),  -- N=48   depth=1
		 4 => (radix => 1, size =>    1),  -- N=60   bypass
		 5 => (radix => 2, size =>    2),  -- N=72   depth=1
		 6 => (radix => 2, size =>    4),  -- N=96   depth=2
		 7 => (radix => 2, size =>    2),  -- N=108  depth=1
		 8 => (radix => 2, size =>    2),  -- N=120  depth=1
		 9 => (radix => 2, size =>    4),  -- N=144  depth=2
		10 => (radix => 2, size =>    2),  -- N=180  depth=1
		11 => (radix => 2, size =>    8),  -- N=192  depth=4
		12 => (radix => 2, size =>    4),  -- N=216  depth=2
		13 => (radix => 2, size =>    4),  -- N=240  depth=2
		14 => (radix => 2, size =>    8),  -- N=288  depth=4
		15 => (radix => 2, size =>    2),  -- N=300  depth=1
		16 => (radix => 2, size =>    4),  -- N=324  depth=2
		17 => (radix => 2, size =>    4),  -- N=360  depth=2
		18 => (radix => 2, size =>   16),  -- N=384  depth=8
		19 => (radix => 2, size =>    8),  -- N=432  depth=4
		20 => (radix => 2, size =>    8),  -- N=480  depth=4
		21 => (radix => 2, size =>    4),  -- N=540  depth=2
		22 => (radix => 2, size =>   16),  -- N=576  depth=8
		23 => (radix => 2, size =>    4),  -- N=600  depth=2
		24 => (radix => 2, size =>    8),  -- N=648  depth=4
		25 => (radix => 2, size =>    8),  -- N=720  depth=4
		26 => (radix => 2, size =>   32),  -- N=768  depth=16
		27 => (radix => 2, size =>   16),  -- N=864  depth=8
		28 => (radix => 2, size =>    4),  -- N=900  depth=2
		29 => (radix => 2, size =>   16),  -- N=960  depth=8
		30 => (radix => 3, size =>   12),  -- N=972  depth=4
		31 => (radix => 2, size =>    8),  -- N=1080 depth=4
		32 => (radix => 2, size =>   32),  -- N=1152 depth=16
		33 => (radix => 2, size =>    8),  -- N=1200 depth=4
		34 => (radix => 2, size =>   16),  -- N=1296 depth=8
		35 => (radix => 2, size =>   16),  -- N=1440 depth=8
		36 => (radix => 2, size =>    4),  -- N=1500 depth=2
		37 => (radix => 2, size =>   64),  -- N=1536 depth=32
		38 => (radix => 3, size =>   12),  -- N=1620 depth=4
		39 => (radix => 2, size =>   32),  -- N=1728 depth=16
		40 => (radix => 2, size =>    8),  -- N=1800 depth=4
		41 => (radix => 2, size =>   32),  -- N=1920 depth=16
		42 => (radix => 3, size =>   24),  -- N=1944 depth=8
		43 => (radix => 2, size =>   16),  -- N=2160 depth=8
		44 => (radix => 2, size =>   64),  -- N=2304 depth=32
		45 => (radix => 2, size =>   16),  -- N=2400 depth=8
		46 => (radix => 2, size =>   32),  -- N=2592 depth=16
		47 => (radix => 3, size =>   12),  -- N=2700 depth=4
		48 => (radix => 2, size =>   32),  -- N=2880 depth=16
		49 => (radix => 3, size =>   36),  -- N=2916 depth=12
		50 => (radix => 2, size =>    8),  -- N=3000 depth=4
		51 => (radix => 2, size =>  128),  -- N=3072 depth=64
		52 => (radix => 3, size =>   24)   -- N=3240 depth=8
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

	function "+" (left, right : t_cmplx_wide) return t_cmplx_wide is
    variable result : t_cmplx_wide;
	begin
		result.re := resize(left.re + right.re, c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, fixed_round);
		result.im := resize(left.im + right.im, c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, fixed_round);
		return result;
	end function;

	function "-" (left, right : t_cmplx_wide) return t_cmplx_wide is
		variable result : t_cmplx_wide;
	begin
		result.re := resize(left.re - right.re, c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, fixed_round);
		result.im := resize(left.im - right.im, c_fxp_int_wide_width-1, -c_fxp_frac_wide_width, fixed_wrap, fixed_round);
		return result;
	end function;

	function "*" (left : t_cmplx_wide; right : t_cmplx_coeff) return t_cmplx_wide is
		variable result 		 : t_cmplx_wide;
		variable mult_result : t_cmplx_wide_mult;
	begin
		mult_result.re := left.re * right.re - left.im * right.im;
		mult_result.im := left.re * right.im + left.im * right.re;
		result.re := resize(mult_result.re, result.re, fixed_wrap, fixed_round);
		result.im := resize(mult_result.im, result.im, fixed_wrap, fixed_round);
		return result;
	end function;

	function "*" (left : t_cmplx; right : t_cmplx_twiddle) return t_cmplx is
		variable result 		 : t_cmplx;
		variable mult_result : t_cmplx_mult;
	begin
		mult_result.re := left.re * right.re - left.im * right.im;
		mult_result.im := left.re * right.im + left.im * right.re;
		result.re := resize(mult_result.re, result.re, fixed_wrap, fixed_round);
		result.im := resize(mult_result.im, result.im, fixed_wrap, fixed_round);
		return result;
	end function;

	function resize(arg : t_cmplx; size_res : t_cmplx_wide) return t_cmplx_wide is
		variable result : t_cmplx_wide;
	begin
		result.re := resize(arg.re, result.re, fixed_wrap, fixed_round);
		result.im := resize(arg.im, result.im, fixed_wrap, fixed_round);
		return result;
	end function;

	function resize(arg : t_cmplx_wide; size_res : t_cmplx) return t_cmplx is
		variable result : t_cmplx;
	begin
		result.re := resize(arg.re, result.re, fixed_wrap, fixed_round);
		result.im := resize(arg.im, result.im, fixed_wrap, fixed_round);
		return result;
	end function;

	function shift_right(arg : t_cmplx_wide; shift_amount : integer) return t_cmplx_wide is
		variable result : t_cmplx_wide;
	begin
		result.re := resize(shift_right(arg.re, shift_amount), result.re, fixed_wrap, fixed_round);
		result.im := resize(shift_right(arg.im, shift_amount), result.im, fixed_wrap, fixed_round);
		return result;
	end function;
 
end package body mr_fft_pkg;