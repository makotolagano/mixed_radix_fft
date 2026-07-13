library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.fixed_pkg.all;
use ieee.fixed_float_types.all;

-- Custom packages
library work;
use work.mr_fft_pkg.all;


entity mr_fft_control is
	generic (
		G_CAPABILITY : natural := 2;
		G_MAX_RADIX : natural := 4
	);
	port (
		i_clk : in  std_logic;
    i_reset : in  std_logic;

		i_phase : in  std_logic_vector(clogb2(G_MAX_RADIX) - 1 downto 0);
		
		-- constants for the current configuration
		o_config_s0 : out std_logic_vector(1 downto 0);
		o_config_s1 : out std_logic;

		-- muxes in front of the preadder
		o_preadder_s0 : out std_logic;
		o_preadder_s1 : out std_logic;

		-- muxes for the input and output of the stage
		o_input_demux_sel : out std_logic_vector(clogb2(G_MAX_RADIX-1) - 1 downto 0);
		o_output_mux_sel : out std_logic_vector(clogb2(G_MAX_RADIX) - 1 downto 0);

		-- FIFOs
		o_fifos_we : out std_logic_vector(0 to G_MAX_RADIX - 2);
		o_fifos_re : out std_logic_vector(0 to G_MAX_RADIX - 2);
		o_fifos_sel : out std_logic_vector(0 to G_MAX_RADIX - 2)
    );
end entity mr_fft_control;

architecture rtl of mr_fft_control is

begin
	
end architecture rtl;