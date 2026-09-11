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

		i_config : in  t_config;

		-- input side phase
		i_phase : in  std_logic_vector(clogb2(G_MAX_RADIX) - 1 downto 0);

		-- preadder mux controls
		o_config_s0 : out std_logic_vector(1 downto 0);
		o_config_s1 : out std_logic;

		-- input demux (which FIFO takes the incoming sample)
		o_input_demux_sel : out std_logic_vector(clogb2(G_MAX_RADIX-1) - 1 downto 0);

		-- delay FIFO write select, one-hot of i_phase. the stage gates it with the input beat.
		o_fifos_we : out std_logic_vector(G_MAX_RADIX - 2 downto 0);
		o_radix_mask : out std_logic_vector(G_MAX_RADIX - 2 downto 0)
    );
end entity mr_fft_control;

architecture rtl of mr_fft_control is
	constant C_FIFO_SEL_ONE : std_logic_vector(G_MAX_RADIX - 2 downto 0) := (std_logic_vector(to_unsigned(1, G_MAX_RADIX - 1)));

	signal fifos_we : std_logic_vector(G_MAX_RADIX - 2 downto 0);
	signal radix_mask : std_logic_vector(G_MAX_RADIX - 2 downto 0);
begin

	-- input demux follows the input phase
	o_input_demux_sel <= i_phase(o_input_demux_sel'length - 1 downto 0);

	-- registered decode (quasi-static)
	PROC_RADIX_MASK: process(i_clk)
	begin
		if rising_edge(i_clk) then
			case i_config.radix is
				when 5 =>
					radix_mask <= std_logic_vector(to_unsigned(2**4-1, G_MAX_RADIX-1));
				when 3 =>
					radix_mask <= std_logic_vector(to_unsigned(2**2-1, G_MAX_RADIX-1));
				when 2 =>
					radix_mask <= std_logic_vector(to_unsigned(1, G_MAX_RADIX-1));
				when others =>
					radix_mask <= (others => '0');
			end case;
		end if;
	end process PROC_RADIX_MASK;

	PROC_FIFO_WE: process(i_phase)
	begin
		fifos_we <= (C_FIFO_SEL_ONE sll to_integer(unsigned(i_phase)));
	end process PROC_FIFO_WE;

	o_fifos_we <= fifos_we;
	o_radix_mask <= radix_mask;

	GEN_PREADDER_MUX_CONTROL_235: if G_CAPABILITY = 2 generate
		-- registered like radix_mask
		PROC_PREADDER_MUX_CONTROL: process(i_clk)
		begin
			if rising_edge(i_clk) then
				if (i_config.radix = 5) then
					o_config_s0 <= "10";
					o_config_s1 <= '1';
				elsif (i_config.radix = 3) then
					o_config_s0 <= "01";
					o_config_s1 <= '0';
				else
					o_config_s0 <= "00";
					o_config_s1 <= '0';
				end if;
			end if;
		end process PROC_PREADDER_MUX_CONTROL;
	end generate GEN_PREADDER_MUX_CONTROL_235;

	GEN_PREADDER_MUX_CONTROL_23: if G_CAPABILITY = 1 generate
		-- registered like radix_mask
		PROC_PREADDER_MUX_CONTROL: process(i_clk)
		begin
			if rising_edge(i_clk) then
				if (i_config.radix = 3) then
					o_config_s0 <= "01";
				else
					o_config_s0 <= "00";
				end if;
			end if;
		end process PROC_PREADDER_MUX_CONTROL;
	end generate GEN_PREADDER_MUX_CONTROL_23;
	
	GEN_PREADDER_MUX_CONTROL_2: if G_CAPABILITY = 0 generate
			o_config_s0 <= (others => '0');
	end generate GEN_PREADDER_MUX_CONTROL_2;
	
end architecture rtl;