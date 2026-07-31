library ieee;
use ieee.std_logic_1164.all;

library work;
use work.mr_fft_pkg.all;
use work.mr_fft_cfg_pkg.all;

-- Synthesis smoke-test wrapper: concrete generics (default: first radix235
-- slot of the mid-bypass pipeline) so ghdl --synth can elaborate the stage.
-- G_SLOT is overridable (xelab -generic_top "G_SLOT=n") to probe any single
-- chain slot in isolation, e.g. when hunting simulator crashes.
entity mr_fft_stage_synth is
  generic (
    G_SLOT : natural := c_num_r2_slots + c_num_r23_slots
  );
  port (
    i_clk   : in  std_logic;
    i_reset : in  std_logic;
    i_config_sel : in std_logic_vector(clogb2(c_num_configs) - 1 downto 0);
    i_sample : in  t_cmplx;
    i_valid  : in  std_logic;
    o_ready  : out std_logic;
    o_sample : out t_cmplx;
    o_valid  : out std_logic;
    i_ready  : in  std_logic
  );
end entity;

architecture rtl of mr_fft_stage_synth is
  constant C_SLOT : natural := G_SLOT;
begin
  dut : entity work.mr_fft_stage
    -- pin the twiddle table into block RAM: under "auto" Vivado's ROM
    -- placement is a heuristic and flip-flops between runs; use "auto"
    -- only for the small radix235 slots (LUT-ROM intended there)
    generic map (G_CAPABILITY => f_slot_capability(C_SLOT), G_CONFIGS => f_slot_configs(C_SLOT),
                 G_TWIDDLE_ROM_STYLE => "block")
    port map (
      i_clk => i_clk, i_reset => i_reset,
      i_config_sel => i_config_sel,
      i_sample => i_sample, i_valid => i_valid, o_ready => o_ready,
      o_sample => o_sample, o_valid => o_valid, i_ready => i_ready);
end architecture;
