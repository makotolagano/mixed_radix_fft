library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.fixed_pkg.all;

library work;
use work.mr_fft_pkg.all;

entity mr_fft_fifo is
	generic (
		G_DEPTH 		 : integer := 1024
	);
	port (
		i_clk 		: in std_logic;
		i_reset 	: in std_logic;

		i_wr_en 	: in std_logic;
		i_wr_sample : in t_cmplx;

		i_rd_en 	: in std_logic;
		o_rd_sample : out t_cmplx;

		o_full 		: out std_logic;
		o_empty 	: out std_logic
	);
end entity mr_fft_fifo;

architecture rtl of mr_fft_fifo is

  constant C_FIFO_DATA_WIDTH : integer := 2*c_fxp_word_width; -- Width of each FIFO data
  signal input_sample : std_logic_vector(C_FIFO_DATA_WIDTH - 1 downto 0);
  signal output_sample : std_logic_vector(C_FIFO_DATA_WIDTH - 1 downto 0);

begin

  input_sample <= to_slv(i_wr_sample.re) & to_slv(i_wr_sample.im);
  o_rd_sample.re <= to_sfixed(output_sample(C_FIFO_DATA_WIDTH - 1 downto C_FIFO_DATA_WIDTH/2), o_rd_sample.re);
  o_rd_sample.im <= to_sfixed(output_sample(C_FIFO_DATA_WIDTH/2 - 1 downto 0), o_rd_sample.im);

	FIFO_INST: entity work.fifo
    generic map (
      G_DATA_WIDTH => C_FIFO_DATA_WIDTH,
      G_DEPTH => G_DEPTH
    )
    port map (
      i_clk 		=> i_clk,
      i_reset 	=> i_reset,

      i_wr_en 	=> i_wr_en,
      i_wr_data => input_sample,

      i_rd_en 	=> i_rd_en,
      o_rd_data => output_sample,

      o_full 		=> o_full,
      o_empty 	=> o_empty
    );

end architecture rtl;
