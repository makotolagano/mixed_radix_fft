library ieee;
use ieee.std_logic_1164.all;

library work;
use work.mr_fft_pkg.all;


-- Materialized mid-bypass pipeline configuration. Every constant here is
-- computed at elaboration time from (c_fft_size_base, c_max_fft_size) in
-- mr_fft_pkg -- the same rule as model/decomposition_configs.py and
-- docs/reconfigurable_fft_midbypass_architecture.md; nothing is hand-written.
--
-- The per-slot config tables are NOT materialized here: instantiate the stage
-- at slot s with G_CONFIGS => f_slot_configs(s) and G_CAPABILITY =>
-- f_slot_capability(s) (both in mr_fft_pkg), as mr_fft_chain does. A deferred
-- package constant of the nested table type crashes xsim 2022.2's debug-
-- database reader (segfault in readAllRTTypes at simulation start).
package mr_fft_cfg_pkg is

	-- number of supported FFT lengths (index space of i_config_sel)
	constant c_num_configs : natural := f_num_configs;
	-- supported FFT lengths, ascending; index = config index
	constant c_fft_sizes : t_nat_arr := f_fft_sizes;

	-- pipeline layout: c_num_r2_slots x radix2 stages, then c_num_r23_slots
	-- x radix23, then c_num_r235_slots x radix235
	constant c_num_r2_slots   : natural := f_num_r2_slots;
	constant c_num_r23_slots  : natural := f_num_r23_slots;
	constant c_num_r235_slots : natural := f_num_r235_slots;
	constant c_num_slots      : natural := f_num_slots;

end package mr_fft_cfg_pkg;
