library ieee;
use ieee.std_logic_1164.all;
use ieee.fixed_pkg.all;
use ieee.fixed_float_types.all;

-- Custom packages
library work;
use work.mr_fft_pkg.all;


entity mr_fft_stage is
	generic (
		G_CAPABILITY : natural := 2
	);
	port (
		i_clk : in  std_logic;
    i_reset : in  std_logic;

		i_sample : in  t_cmplx;

		i_s0 : in  std_logic_vector(1 downto 0);
		i_s1 : in  std_logic;

		o_sample : out t_cmplx
	);
end entity mr_fft_stage;

architecture rtl of mr_fft_stage is

  function get_max_radix(g_capability : natural) return natural is
  begin
    if g_capability = 2 then
      return 5;
    elsif g_capability = 1 then
      return 3;
    elsif g_capability = 0 then
      return 2;
    else
      return 0;
    end if;
  end function get_max_radix;

  constant C_MAX_RADIX       : natural := get_max_radix(G_CAPABILITY);
  constant C_NUM_FIFOS       : natural := C_MAX_RADIX - 1; -- Number of FIFOs needed for the given capability
  constant C_FIFO_DATA_WIDTH : natural := 36; -- Width of each FIFO data
  constant C_FIFO_DEPTH      : natural := 1024;    -- Depth of each FIFO

  signal input_sample : t_cmplx;

  type t_fifo_data_array is array (0 to C_NUM_FIFOS - 1) of t_cmplx;
  signal fifos_data_in  : t_fifo_data_array;
  signal fifos_data_out : t_fifo_data_array;
  signal fifos_we       : std_logic_vector(0 to C_NUM_FIFOS - 1);
  signal fifos_re       : std_logic_vector(0 to C_NUM_FIFOS - 1);
  signal fifos_sel      : std_logic_vector(0 to C_NUM_FIFOS - 1);

  type t_preadder_signals_array is array (0 to C_MAX_RADIX - 1) of t_cmplx;
  signal preadder_inputs  : t_preadder_signals_array;
  signal preadder_outputs : t_preadder_signals_array;
  signal preadder_input2_muxed : t_cmplx;

  signal s0 : std_logic_vector(1 downto 0);
  signal s1 : std_logic;
  signal preadder_s0 : std_logic;
  signal preadder_s1 : std_logic;

  type t_input_demux_array is array (0 to C_NUM_FIFOS - 1) of t_cmplx;
  signal input_demux_out : t_input_demux_array;
  signal input_demux_sel : std_logic_vector(clogb2(C_NUM_FIFOS) - 1 downto 0);

  type t_output_mux_array is array (0 to C_MAX_RADIX - 1) of t_cmplx;
  signal output_mux_in : t_output_mux_array;
  signal output_mux_out : t_cmplx;
  signal output_mux_sel : std_logic_vector(clogb2(C_MAX_RADIX-1) - 1 downto 0);

begin

  input_sample <= i_sample;

	GEN_FIFOS: for i in 0 to C_NUM_FIFOS - 1 generate
    FIFO_INST: entity work.mr_fft_fifo
      generic map (
        G_DEPTH      => C_FIFO_DEPTH
      )
      port map (
        i_clk     => i_clk,
        i_reset   => i_reset,
        i_wr_en   => fifos_we(i),
        i_wr_sample => fifos_data_in(i),
        i_rd_en   => fifos_re(i),
        o_rd_sample => fifos_data_out(i)
      );

      
  end generate GEN_FIFOS;

  GEN_INPUTS_TO_FIFOS: for i in 0 to C_NUM_FIFOS - 1 generate
    fifos_data_in(i) <= preadder_outputs(i+1) when fifos_sel(i) = '1' else input_demux_out(i);
  end generate GEN_INPUTS_TO_FIFOS;

  -- Input demux logic

  GEN_INPUT_DEMUX_235: if G_CAPABILITY = 2 generate
    PROC_INPUT_DEMUX: process(input_sample, input_demux_sel)
    begin
      input_demux_out <= (others => (others => "0")); -- Default assignment
      case input_demux_sel is
        when "00" =>
          input_demux_out(0) <= input_sample;
        when "01" =>
          input_demux_out(1) <= input_sample;
        when "10" =>
          input_demux_out(2) <= input_sample;
        when "11" =>
          input_demux_out(3) <= input_sample;
      end case;
    end process PROC_INPUT_DEMUX;
  end generate GEN_INPUT_DEMUX_235;

  GEN_INPUT_DEMUX_23: if G_CAPABILITY = 1 generate
    PROC_INPUT_DEMUX: process(input_sample, input_demux_sel)
    begin
      input_demux_out <= (others => (others => "0")); -- Default assignment
      case input_demux_sel is
        when "0" =>
          input_demux_out(0) <= input_sample;
        when "1" =>
          input_demux_out(1) <= input_sample;
      end case;
    end process PROC_INPUT_DEMUX;
  end generate GEN_INPUT_DEMUX_23;

  GEN_INPUT_DEMUX_2: if G_CAPABILITY = 0 generate
    input_demux_out(0) <= input_sample;
  end generate GEN_INPUT_DEMUX_2;
  
  GEN_PREADDER_INPUTS_235: if G_CAPABILITY = 2 generate
    preadder_s0 <= not s0(1);
    preadder_s1 <= not s1;

    preadder_inputs(0) <= fifos_data_out(0);
    preadder_inputs(1) <= input_sample when preadder_s0 = '1' else fifos_data_out(1);
    
    preadder_input2_muxed <= input_sample when preadder_s1 = '1' else fifos_data_out(2);
    preadder_inputs(2) <= preadder_input2_muxed when not(preadder_s0 = '1' and preadder_s1 = '1') else (others => (others => '0'));
    preadder_inputs(3) <= fifos_data_out(3) when preadder_s1 = '0' else (others => (others => '0'));
    preadder_inputs(4) <= input_sample when preadder_s1 = '0' else (others => (others => '0'));
  end generate GEN_PREADDER_INPUTS_235;

  GEN_PREADDER_INPUTS_23: if G_CAPABILITY = 1 generate
    preadder_s0 <= not s0(0);

    preadder_inputs(0) <= fifos_data_out(0);
    preadder_inputs(1) <= input_sample when preadder_s0 = '1' else fifos_data_out(1);
    preadder_inputs(2) <= input_sample when preadder_s0 = '0' else (others => (others => '0'));
  end generate GEN_PREADDER_INPUTS_23;

  GEN_PREADDER_INPUTS_2: if G_CAPABILITY = 0 generate
    preadder_inputs(0) <= fifos_data_out(0);
    preadder_inputs(1) <= input_sample;
  end generate GEN_PREADDER_INPUTS_2;

  PREADDER_INST: entity work.mr_fft_preadder
    generic map (
      G_CAPABILITY => G_CAPABILITY
    )
    port map (
      i_clk => i_clk,

      i_x0 => preadder_inputs(0),
      i_x1 => preadder_inputs(1),
      i_x2 => preadder_inputs(2),
      i_x3 => preadder_inputs(3),
      i_x4 => preadder_inputs(4),

      i_s0 => s0,
      i_s1 => s1,

      o_X0 => preadder_outputs(0),
      o_X1 => preadder_outputs(1),
      o_X2 => preadder_outputs(2),
      o_X3 => preadder_outputs(3),
      o_X4 => preadder_outputs(4)
    );

  GEN_MUX_OUTPUT_235: if G_CAPABILITY = 2 generate
    PROC_OUTPUT: process(preadder_outputs(0), fifos_data_out, output_mux_sel)
    begin
      output_mux_out <= (others => (others => '0')); -- Default assignment
      case output_mux_sel is
        when "000" =>
          output_mux_out <= fifos_data_out(0);
        when "001" =>
          output_mux_out <= fifos_data_out(1);
        when "010" =>
          output_mux_out <= fifos_data_out(2);
        when "011" =>
          output_mux_out <= fifos_data_out(3);
        when "100" =>
          output_mux_out <= preadder_outputs(0);
      end case;
    end process PROC_OUTPUT;
	end generate GEN_MUX_OUTPUT_235;

  GEN_MUX_OUTPUT_23: if G_CAPABILITY = 1 generate
    PROC_OUTPUT: process(preadder_outputs(0), fifos_data_out, output_mux_sel)
    begin
      output_mux_out <= (others => (others => '0')); -- Default assignment
      case output_mux_sel is
        when "00" =>
          output_mux_out <= fifos_data_out(0);
        when "01" =>
          output_mux_out <= fifos_data_out(1);
        when "10" =>
          output_mux_out <= preadder_outputs(0);
      end case;
    end process PROC_OUTPUT;
  end generate GEN_MUX_OUTPUT_23;

  GEN_MUX_OUTPUT_2: if G_CAPABILITY = 0 generate
    output_mux_out <= fifos_data_out(0) when output_mux_sel = "0" else preadder_outputs(0);
  end generate GEN_MUX_OUTPUT_2;

  o_sample <= output_mux_out;

  
end architecture rtl;