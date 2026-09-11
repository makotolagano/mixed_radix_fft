library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.fixed_pkg.all;

library work;
use work.mr_fft_pkg.all;
use work.mr_fft_cfg_pkg.all;

-- Chain of c_num_slots stages connected through ready/valid.
--
--   [ radix2 x c_num_r2_slots ][ radix23 x c_num_r23_slots ][ radix235 x c_num_r235_slots ]
--
-- i_config_sel selects the FFT length, every stage looks up its own
-- (radix, delay, bypass) from its G_CONFIGS table. a bypassed stage still
-- registers the stream through its output skid, so there are no long
-- combinational paths across slots.
--
-- reconfigure only when drained: feed whole frames, wait for all outputs,
-- then change i_config_sel. no reset needed.
--
-- every stage shifts its output right by ceil(log2(radix)), so the chain
-- output is X * 2**(-total_shift). the final scaler in the top fixes the gain.
entity mr_fft_chain is
	generic (
		G_PIPELINE : boolean := true
	);
	port (
		i_clk   : in  std_logic;
		i_reset : in  std_logic;

		-- FFT length select, change only when drained
		i_config_sel : in std_logic_vector(clogb2(c_num_configs) - 1 downto 0);

		-- input stream handshake
		i_sample : in  t_cmplx;
		i_valid  : in  std_logic;
		o_ready  : out std_logic;

		-- output stream handshake
		o_sample : out t_cmplx;
		o_valid  : out std_logic;
		i_ready  : in  std_logic
	);
end entity mr_fft_chain;

architecture rtl of mr_fft_chain is

	-- inter-slot links: link 0 = chain input, link c_num_slots = chain output
	type t_link_data is array (0 to c_num_slots) of t_cmplx;
	signal link_data  : t_link_data;
	signal link_valid : std_logic_vector(0 to c_num_slots);
	signal link_ready : std_logic_vector(0 to c_num_slots);

begin

	-- chain boundaries
	link_data(0)  <= i_sample;
	link_valid(0) <= i_valid;
	o_ready       <= link_ready(0);

	o_sample                 <= link_data(c_num_slots);
	o_valid                  <= link_valid(c_num_slots);
	link_ready(c_num_slots)  <= i_ready;

	GEN_SLOTS: for s in 0 to c_num_slots - 1 generate
		constant C_CAP  : natural := f_slot_capability(s);
		constant C_CFGS : t_config_arr(0 to c_num_configs - 1) := f_slot_configs(s);
		-- big twiddle tables go to block RAM, the small radix235 ones stay in LUTs
		constant C_ROM_BLOCK : boolean := C_CAP /= 2;
		constant C_FIFO_RAM_STYLE : string := f_ram_style(get_delay_cnt(C_CFGS));
	begin

		GEN_STAGE_BRAM_ROM: if C_ROM_BLOCK generate
			STAGE_INST: entity work.mr_fft_stage
				generic map (
					G_CAPABILITY        => C_CAP,
					G_CONFIGS           => C_CFGS,
					G_PIPELINE          => G_PIPELINE,
					G_TWIDDLE_ROM_STYLE => "block",
					G_FIFO_RAM_STYLE    => C_FIFO_RAM_STYLE
				)
				port map (
					i_clk    => i_clk,
					i_reset  => i_reset,
					i_config_sel => i_config_sel,
					i_sample => link_data(s),
					i_valid  => link_valid(s),
					o_ready  => link_ready(s),
					o_sample => link_data(s + 1),
					o_valid  => link_valid(s + 1),
					i_ready  => link_ready(s + 1)
				);
		end generate GEN_STAGE_BRAM_ROM;

		GEN_STAGE_LUT_ROM: if not C_ROM_BLOCK generate
			STAGE_INST: entity work.mr_fft_stage
				generic map (
					G_CAPABILITY        => C_CAP,
					G_CONFIGS           => C_CFGS,
					G_PIPELINE          => G_PIPELINE,
					G_TWIDDLE_ROM_STYLE => "auto",
					G_FIFO_RAM_STYLE    => C_FIFO_RAM_STYLE
				)
				port map (
					i_clk    => i_clk,
					i_reset  => i_reset,
					i_config_sel => i_config_sel,
					i_sample => link_data(s),
					i_valid  => link_valid(s),
					o_ready  => link_ready(s),
					o_sample => link_data(s + 1),
					o_valid  => link_valid(s + 1),
					i_ready  => link_ready(s + 1)
				);
		end generate GEN_STAGE_LUT_ROM;

	end generate GEN_SLOTS;

end architecture rtl;
