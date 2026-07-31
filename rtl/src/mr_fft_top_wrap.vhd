library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.fixed_pkg.all;

library work;
-- ---------------------------------------------------------------------------
-- Mixed-radix FFT top: the mid-bypass chain plus an AXI4-Lite configuration
-- interface. Single clock domain (AXI and stream share i_clk).
--
-- Register map (32-bit registers, byte addresses):
--   0x00  ID           RO  x"0FF70100" (FFT core, v1.0)
--   0x04  CTRL         WO  bit0 COMMIT (self-clearing), bit1 CLR_FERR
--   0x08  STATUS       RO  bit0 BUSY (commit pending), bit1 IDLE (in_flight=0),
--                          bit3 FRAMING_ERR (sticky, clear via CTRL.CLR_FERR)
--   0x0C  CONFIG_SEL   RW  SHADOW config select (index into c_fft_sizes)
--   0x10  CONFIG_ACTIVE RO active (committed) config select
--   0x14  IN_FLIGHT    RO  input beats minus output beats inside the chain
--   0x18  FFT_SIZE     RO  N of the ACTIVE config (from c_fft_sizes)
--
-- Configuration protocol (encodes the drained-switch contract in hardware):
-- writes to CONFIG_SEL land in a shadow register only. Writing CTRL.COMMIT
-- gates the input stream (o_ready forced low), waits until the chain is
-- drained (IN_FLIGHT = 0 -- exact, by the chain's beat conservation), then
-- transfers shadow -> active in one cycle and reopens the input. Software:
-- write CONFIG_SEL, write COMMIT, poll STATUS.BUSY = 0 (or just keep
-- streaming: the gate handles the boundary). NOTE: commit only after WHOLE
-- frames have been offered -- a partial frame parks samples in the delay
-- FIFOs and IN_FLIGHT never reaches zero (STATUS makes this visible).
--
-- Frame marker (AXI-Stream tlast semantics): o_last is REGENERATED from a
-- mod-N output beat counter (N of the ACTIVE config), never forwarded from
-- the input -- exact because a commit only lands with the chain drained on a
-- frame boundary. i_last cannot steer the datapath (the delay FIFOs are
-- configured for N), so it is only checked: an accepted beat where i_last
-- disagrees with the input beat counter sets the sticky STATUS.FRAMING_ERR.
-- ---------------------------------------------------------------------------
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
