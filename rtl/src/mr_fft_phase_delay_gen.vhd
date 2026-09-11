library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.fixed_pkg.all;
use ieee.fixed_float_types.all;

-- Custom packages
library work;
use work.mr_fft_pkg.all;


entity mr_fft_phase_delay_gen is
	generic (
		G_CAPABILITY : natural := 2;
		G_DELAY_CNT : natural := 16;
		G_MAX_RADIX : natural := 4
	);
	port (
		i_clk : in  std_logic;
    i_reset : in  std_logic;

    i_config_radix : in  std_logic_vector(clogb2(G_MAX_RADIX) - 1 downto 0);
		i_config_delay : in  std_logic_vector(clogb2(G_DELAY_CNT) - 1 downto 0);
		i_en : in  std_logic;

    o_phase : out std_logic_vector(clogb2(G_MAX_RADIX) - 1 downto 0);
		o_delay_cnt : out std_logic_vector(clogb2(G_DELAY_CNT) - 1 downto 0)
    );
end entity mr_fft_phase_delay_gen;

architecture rtl of mr_fft_phase_delay_gen is
	signal phase_cnt : unsigned(clogb2(G_MAX_RADIX) - 1 downto 0) := (others => '0');
	signal delay_cnt : unsigned(clogb2(G_DELAY_CNT) - 1 downto 0) := (others => '0');

	-- registered -1 values so the wrap compares start from flops
	signal radix_m1 : unsigned(clogb2(G_MAX_RADIX) - 1 downto 0);
	signal delay_m1 : unsigned(clogb2(G_DELAY_CNT) - 1 downto 0);
begin

	PROC_CFG_M1: process(i_clk)
	begin
		if rising_edge(i_clk) then
			radix_m1 <= unsigned(i_config_radix) - 1;
			delay_m1 <= unsigned(i_config_delay) - 1;
		end if;
	end process PROC_CFG_M1;

	PROC_PHASE_CNT: process(i_clk)
	begin
		if rising_edge(i_clk) then
			if (i_reset = '1') then
				phase_cnt <= (others => '0');
			else
				if (i_en = '1' and (delay_cnt = delay_m1)) then
					phase_cnt <= phase_cnt + 1;
					if (phase_cnt = radix_m1) then
						phase_cnt <= (others => '0');
					end if;
				end if;
			end if;
		end if;
	end process PROC_PHASE_CNT;

	PROC_DELAY_CNT: process(i_clk)
	begin
		if rising_edge(i_clk) then
			if (i_reset = '1') then
				delay_cnt <= (others => '0');
			else
				if (i_en = '1') then
					delay_cnt <= delay_cnt + 1;
					if (delay_cnt = delay_m1) then
						delay_cnt <= (others => '0');
					end if;
				end if;
			end if;
		end if;
	end process PROC_DELAY_CNT;

	o_phase <= std_logic_vector(phase_cnt);
	o_delay_cnt <= std_logic_vector(delay_cnt);
	
end architecture rtl;