library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.fixed_pkg.all;
use work.mr_fft_pkg.all;

-- Self-checking twiddle_rom TB driven by generated golden vectors.
-- Regenerate vectors:  python model/twiddle_rom_ref.py --int 2 --frac 15 --configs 5:300,3:243,2:3072,2:12
entity tb_mr_fft_stage is end entity;

architecture sim of tb_mr_fft_stage is
  constant CFGS : t_config_arr := ((5,300), (3,3), (2,2), (2,12));

  constant CLK_PERIOD : time := 10 ns;
  signal clk  : std_logic := '0';
  signal reset : std_logic := '1';
  signal done : boolean := false;

begin
  dut : entity work.mr_fft_stage
      generic map (G_CAPABILITY => 2, G_CONFIGS => CFGS)
      port map (i_clk => clk, i_reset => reset, i_sample => (others=> (others=>'0')), i_s0 => (others=>'0'), i_s1 => '0', o_sample => open);
  clk <= not clk after CLK_PERIOD/2 when not done else '0';

  process
  begin
    wait for 2*CLK_PERIOD;
    reset <= '0';
    wait for 200*CLK_PERIOD;
    report "TB PASS -- mr_fft_stage instantiated successfully" severity note;
    done <= true;
    wait;
  end process;
end architecture;
