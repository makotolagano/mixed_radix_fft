library ieee;
use ieee.std_logic_1164.all;

library work;
use work.mr_fft_pkg.all;


-- Pipeline configuration, all computed at elaboration from
-- (c_fft_size_base, c_max_fft_size) in mr_fft_pkg.
--
-- per-slot config tables are not stored here, instantiate the stage with
-- G_CONFIGS => f_slot_configs(s). a deferred constant of that type crashes xsim 2022.2.
package mr_fft_cfg_pkg is

	-- number of supported FFT lengths
	constant c_num_configs : natural := f_num_configs;
	-- supported FFT lengths, ascending; index = config index
	constant c_fft_sizes : t_nat_arr := f_fft_sizes;
	-- final scaler codes, index aligned
	constant c_fft_scales : t_nat_arr := f_fft_scales;

	-- slot layout: radix2 slots, then radix23, then radix235
	constant c_num_r2_slots   : natural := f_num_r2_slots;
	constant c_num_r23_slots  : natural := f_num_r23_slots;
	constant c_num_r235_slots : natural := f_num_r235_slots;
	constant c_num_slots      : natural := f_num_slots;

end package mr_fft_cfg_pkg;
