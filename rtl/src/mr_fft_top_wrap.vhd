library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.fixed_pkg.all;

library work;
-- Active-low reset wrapper around mr_fft_top for the block design.
-- register map and protocol: see mr_fft_top.
entity mr_fft_top_wrap is
	port (
		i_clk   : in  std_logic;
		i_reset_n : in  std_logic;

		-- AXI4-Lite slave (configuration)
		s_axi_awaddr  : in  std_logic_vector(7 downto 0);
		s_axi_awvalid : in  std_logic;
		s_axi_awready : out std_logic;
		s_axi_wdata   : in  std_logic_vector(31 downto 0);
		s_axi_wstrb   : in  std_logic_vector(3 downto 0);
		s_axi_wvalid  : in  std_logic;
		s_axi_wready  : out std_logic;
		s_axi_bresp   : out std_logic_vector(1 downto 0);
		s_axi_bvalid  : out std_logic;
		s_axi_bready  : in  std_logic;
		s_axi_araddr  : in  std_logic_vector(7 downto 0);
		s_axi_arvalid : in  std_logic;
		s_axi_arready : out std_logic;
		s_axi_rdata   : out std_logic_vector(31 downto 0);
		s_axi_rresp   : out std_logic_vector(1 downto 0);
		s_axi_rvalid  : out std_logic;
		s_axi_rready  : in  std_logic;

		-- input stream handshake
		i_sample : in  std_logic_vector(63 downto 0);
		i_valid  : in  std_logic;
		i_last   : in  std_logic;
		o_ready  : out std_logic;

		-- output stream handshake
		o_sample : out std_logic_vector(63 downto 0);
		o_valid  : out std_logic;
		o_last   : out std_logic;
		i_ready  : in  std_logic
	);
end entity mr_fft_top_wrap;

architecture rtl of mr_fft_top_wrap is
    signal reset : std_logic;

begin
  reset <= not i_reset_n;

  mr_fft_top_inst : entity work.mr_fft_top
		generic map (G_PIPELINE => true)
		port map (
			i_clk => i_clk, i_reset => reset,
			s_axi_awaddr => s_axi_awaddr, s_axi_awvalid => s_axi_awvalid, s_axi_awready => s_axi_awready,
			s_axi_wdata => s_axi_wdata, s_axi_wstrb => s_axi_wstrb, s_axi_wvalid => s_axi_wvalid,
			s_axi_wready => s_axi_wready, s_axi_bresp => s_axi_bresp, s_axi_bvalid => s_axi_bvalid,
			s_axi_bready => s_axi_bready,
			s_axi_araddr => s_axi_araddr, s_axi_arvalid => s_axi_arvalid, s_axi_arready => s_axi_arready,
			s_axi_rdata => s_axi_rdata, s_axi_rresp => s_axi_rresp, s_axi_rvalid => s_axi_rvalid,
			s_axi_rready => s_axi_rready,
			i_sample => i_sample, i_valid => i_valid, i_last => i_last,
			o_ready => o_ready,
			o_sample => o_sample, o_valid => o_valid, o_last => o_last,
			i_ready => i_ready);

end architecture rtl;
