library ieee;
use ieee.std_logic_1164.all;
use ieee.fixed_pkg.all;
use ieee.fixed_float_types.all;

-- Custom packages
library work;
use work.mr_fft_pkg.all;


entity mr_fft_control is
	generic (
		G_CAPABILITY : natural := 2
	);
	port (
		i_clk : in  std_logic;
    i_reset : in  std_logic;

    i_config : in  std_logic_vector(1 downto 0);
    o_phase : out std_logic_vector(1 downto 0)
    );
end entity mr_fft_control;

architecture rtl of mr_fft_control is

begin
	
end architecture rtl;