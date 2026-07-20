library ieee;
use ieee.std_logic_1164.all;

library work;
use work.mr_fft_pkg.all;


-- Materialized mid-bypass pipeline configuration. Every constant here is
-- computed at elaboration time from (c_fft_size_base, c_max_fft_size) in
-- mr_fft_pkg -- the same rule as model/decomposition_configs.py and
-- docs/reconfigurable_fft_midbypass_architecture.md; nothing is hand-written.
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

	type t_slot_config_arr is array (natural range <>) of t_config_arr(0 to c_num_configs - 1);

	-- per-slot config tables: c_slot_cfgs(s)(c) = (radix, size) of config c
	-- at pipeline slot s -- (1, 1) where the slot is bypassed. Instantiate
	-- the stage at slot s with G_CONFIGS => c_slot_cfgs(s) and
	-- G_CAPABILITY => f_slot_capability(s).
	constant c_slot_cfgs : t_slot_config_arr;

end package mr_fft_cfg_pkg;


package body mr_fft_cfg_pkg is

	function f_all_slot_cfgs return t_slot_config_arr is
		variable cfgs : t_slot_config_arr(0 to c_num_slots - 1);
	begin
		for s in cfgs'range loop
			cfgs(s) := f_slot_configs(s);
		end loop;
		return cfgs;
	end function f_all_slot_cfgs;

	constant c_slot_cfgs : t_slot_config_arr := f_all_slot_cfgs;

end package body mr_fft_cfg_pkg;
