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

		i_phase : in  std_logic_vector(clogb2(G_MAX_RADIX) - 1 downto 0);
		
		-- constants for the current configuration
		o_config_s0 : out std_logic_vector(1 downto 0);
		o_config_s1 : out std_logic;

		-- muxes for the input and output of the stage
		o_input_demux_sel : out std_logic_vector(clogb2(G_MAX_RADIX-1) - 1 downto 0);
		o_output_mux_sel : out std_logic_vector(clogb2(G_MAX_RADIX) - 1 downto 0);

		-- FIFOs
		o_fifos_we : out std_logic_vector(G_MAX_RADIX - 2 downto 0);
		o_fifos_re : out std_logic_vector(G_MAX_RADIX - 2 downto 0);
		o_fifos_sel : out std_logic_vector(G_MAX_RADIX - 2 downto 0)
    );
end entity mr_fft_control;

architecture rtl of mr_fft_control is
	constant C_FIFO_SEL_ONE : std_logic_vector(G_MAX_RADIX - 2 downto 0) := (std_logic_vector(to_unsigned(1, G_MAX_RADIX - 1)));

	signal fifos_we : std_logic_vector(G_MAX_RADIX - 2 downto 0);
begin

	-- Input demux and output mux signals based on the phase
	o_input_demux_sel <= i_phase(o_input_demux_sel'length - 1 downto 0);

	PROC_FIFO_SEL: process(i_config, i_phase)
	begin
		if (i_config.radix-1 = to_integer(unsigned(i_phase))) then
			o_fifos_sel <= (others => '1');
		else
			o_fifos_sel <= (others => '0');
		end if;
	end process PROC_FIFO_SEL;

	PROC_FIFO_WE: process(i_config, i_phase)
	begin
		if (i_config.radix-1 = to_integer(unsigned(i_phase))) then
			fifos_we <= (others => '1');
		else
			fifos_we <= (C_FIFO_SEL_ONE sll to_integer(unsigned(i_phase)));
		end if;
	end process PROC_FIFO_WE;

	o_fifos_we <= fifos_we;
	o_fifos_re <= fifos_we;

	GEN_OUTPUT_MUX_SEL_235: if G_CAPABILITY = 2 generate
		-- Output the output mux signals based on the phase
		PROC_OUTPUT_MUX_SEL: process(i_config, i_phase)
		begin
			if (i_config.radix = 5) then
				o_output_mux_sel <= i_phase;
			elsif (i_config.radix = 3) then
				o_output_mux_sel <= i_phase(i_phase'length - 1) & '0' & i_phase(0);
			else
				o_output_mux_sel <= i_phase(i_phase'length - 1) & "00";
			end if;
		end process PROC_OUTPUT_MUX_SEL;
	end generate GEN_OUTPUT_MUX_SEL_235;

	GEN_OUTPUT_MUX_SEL_23: if G_CAPABILITY = 1 generate
		-- Output the output mux signals based on the phase
		PROC_OUTPUT_MUX_SEL: process(i_config, i_phase)
		begin
			if (i_config.radix = 3) then
				o_output_mux_sel <= i_phase;
			else
				o_output_mux_sel <= i_phase(i_phase'length - 1) & '0';
			end if;
		end process PROC_OUTPUT_MUX_SEL;
	end generate GEN_OUTPUT_MUX_SEL_23;

	GEN_OUTPUT_MUX_SEL_2: if G_CAPABILITY = 0 generate
		-- Output the output mux signals based on the phase
		PROC_OUTPUT_MUX_SEL: process(i_phase)
		begin
			o_output_mux_sel <= i_phase;
		end process PROC_OUTPUT_MUX_SEL;
	end generate GEN_OUTPUT_MUX_SEL_2;

	GEN_PREADDER_MUX_CONTROL_235: if G_CAPABILITY = 2 generate
		-- Output the preadder mux signals based on the configuration
		PROC_PREADDER_MUX_CONTROL: process(i_config)
		begin
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
		end process PROC_PREADDER_MUX_CONTROL;
	end generate GEN_PREADDER_MUX_CONTROL_235;
	
	GEN_PREADDER_MUX_CONTROL_23: if G_CAPABILITY = 1 generate
		-- Output the preadder mux signals based on the configuration
		PROC_PREADDER_MUX_CONTROL: process(i_config)
		begin
			if (i_config.radix = 3) then
				o_config_s0 <= "1";
			else
				o_config_s0 <= "0";
			end if;
		end process PROC_PREADDER_MUX_CONTROL;
	end generate GEN_PREADDER_MUX_CONTROL_23;
	
	GEN_PREADDER_MUX_CONTROL_2: if G_CAPABILITY = 0 generate
			o_config_s0 <= (others => '0');
	end generate GEN_PREADDER_MUX_CONTROL_2;
	
end architecture rtl;