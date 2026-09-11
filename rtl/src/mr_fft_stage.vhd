library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.fixed_pkg.all;
use ieee.fixed_float_types.all;

-- Custom packages
library work;
use work.mr_fft_pkg.all;


entity mr_fft_stage is
	generic (
		G_CAPABILITY : natural := 2;
    G_CONFIGS : t_config_arr;
    -- pipeline the preadder; the stage runs a valid pipeline of the same depth.
    -- ignored for capability 0, the radix-2 butterfly stays combinational.
    G_PIPELINE : boolean := true;
    -- rom_style of the twiddle table, "block" for the big tables
    G_TWIDDLE_ROM_STYLE : string := "auto";
    G_FIFO_RAM_STYLE    : string := "auto"
	);
	port (
		i_clk : in  std_logic;
    i_reset : in  std_logic;

    -- index into G_CONFIGS. radix, delay, twiddle block and bypass are all looked up from it.
    i_config_sel : in std_logic_vector(clogb2(G_CONFIGS'length) - 1 downto 0) := (others => '0');

    -- no reset needed to reconfigure: feed whole frames, drain all outputs, then change
    -- i_config_sel. o_ready is held low for a few cycles after a switch.

    -- input and output sides are decoupled. delay FIFOs hold input samples only,
    -- preadder results land in result FIFOs and are drained when downstream is ready.
    -- the input side never looks at i_ready.

    -- input stream handshake
		i_sample : in  t_cmplx;
    i_valid  : in  std_logic := '1';
    o_ready  : out std_logic;

    -- output stream handshake
		o_sample : out t_cmplx;
    o_valid  : out std_logic;
    i_ready  : in  std_logic := '1'
	);
end entity mr_fft_stage;

architecture rtl of mr_fft_stage is

  constant C_MAX_RADIX       : natural := get_max_radix(G_CAPABILITY);
  constant C_NUM_FIFOS       : natural := C_MAX_RADIX - 1; -- Number of FIFOs needed for the given capability
  constant C_FIFO_DATA_WIDTH : natural := c_fxp_word_width; -- Width of each FIFO data
  constant C_DELAY_CNT       : natural := get_delay_cnt(G_CONFIGS);    -- largest delay (size / radix)
  constant C_NUM_CONFIGS     : natural := G_CONFIGS'length; -- Number of configurations
  constant C_FIFO_DEPTH      : natural := C_DELAY_CNT;    -- Depth of each delay FIFO
  -- capability 0 is never pipelined
  constant C_PIPELINE        : boolean := G_PIPELINE and G_CAPABILITY /= 0;
  -- preadder latency, the valid pipeline matches it
  constant C_PRE_LAT         : natural := preadder_latency(G_CAPABILITY, C_PIPELINE);
  -- one result block plus the pipeline overlap, so streaming never stalls at the joint phase
  constant C_RESULT_DEPTH    : natural := C_DELAY_CNT + C_PRE_LAT + 2;

  -- last stage (delay 1) only rotates by W^0: no ROM, no multiplier
  constant C_NEED_ROT : boolean := C_DELAY_CNT > 1;

  -- block RAM tables use the BRAM output register, one more cycle of latency
  constant C_TW_OUT_REG : boolean := G_TWIDDLE_ROM_STYLE = "block";
  -- twiddle latency after a rot beat: folds + address reg + read + reconstruct
  constant C_TW_LAT     : natural := 5 + boolean'pos(C_TW_OUT_REG);
  -- rotator follows G_PIPELINE directly, capability 0 has a rotator too
  constant C_ROT_LAT    : natural := rotator_latency(G_PIPELINE);

  -- rotator: twiddle exponent range
  constant C_MAX_SIZE : natural := get_max_size(G_CONFIGS);
  constant C_TW_K_W   : natural := clogb2(C_MAX_SIZE);
  -- output skid: the credit counts everything in flight through the rotator
  constant C_SKID_DEPTH : natural := C_TW_LAT + C_ROT_LAT + 4;

  signal input_sample : t_cmplx;

  -- per-config values, all looked up from i_config_sel.
  -- delay uses modulo encoding: delay = 2**width reads as 0, the -1 compares still work.
  constant C_NUM_CFGS : natural := G_CONFIGS'length;

  type t_delay_tbl is array (0 to C_NUM_CFGS - 1) of
    std_logic_vector(clogb2(C_DELAY_CNT) - 1 downto 0);

  function f_delay_tbl return t_delay_tbl is
    variable t : t_delay_tbl;
    variable c : t_config;
  begin
    for i in 0 to C_NUM_CFGS - 1 loop
      c := G_CONFIGS(G_CONFIGS'low + i);
      t(i) := std_logic_vector(to_unsigned(
          (c.size / c.radix) mod 2**(clogb2(C_DELAY_CNT)), clogb2(C_DELAY_CNT)));
    end loop;
    return t;
  end function f_delay_tbl;

  constant C_DELAY_TBL : t_delay_tbl := f_delay_tbl;
  constant C_FIFO_RAM_STYLE : string := f_ram_style(C_DELAY_CNT);

  signal config_sel : std_logic_vector(clogb2(G_CONFIGS'length) - 1 downto 0);
  signal cfg_idx      : natural range 0 to C_NUM_CFGS - 1;
  -- init so radix is never 0 at time zero
  signal config       : t_config := G_CONFIGS(G_CONFIGS'low);
  signal bypass       : std_logic;     -- selected entry is a bypass
  signal config_delay : std_logic_vector(clogb2(C_DELAY_CNT) - 1 downto 0);
  -- delay - 1, registered so the compares start from a flop
  signal delay_last_r : unsigned(clogb2(C_DELAY_CNT) - 1 downto 0);
  -- radix - 1, same reason
  signal radix_m1     : unsigned(clogb2(C_MAX_RADIX) - 1 downto 0);

  -- hold o_ready low a few cycles after i_config_sel changes
  constant C_CFG_SETTLE : natural := 3;
  signal cfg_settled : std_logic_vector(C_CFG_SETTLE - 1 downto 0);
  signal cfg_ok      : std_logic;
  signal o_ready_pre : std_logic;   -- ready before the settle gate

  -- delay FIFOs, input samples only (capability 0 also writes X1 back)
  type t_fifo_data_array is array (0 to C_NUM_FIFOS - 1) of t_cmplx;
  signal fifos_data_in  : t_fifo_data_array;
  signal fifos_data_out : t_fifo_data_array;   -- show-ahead FIFO heads
  signal fifos_we       : std_logic_vector(C_NUM_FIFOS - 1 downto 0);
  signal fifos_rd_valid : std_logic_vector(C_NUM_FIFOS - 1 downto 0);
  signal fifos_full     : std_logic_vector(C_NUM_FIFOS - 1 downto 0);

  -- result FIFOs, one per butterfly output
  type t_result_data_array is array (0 to C_MAX_RADIX - 1) of t_cmplx;
  signal results_data_out : t_result_data_array;   -- show-ahead heads
  signal results_we       : std_logic_vector(C_MAX_RADIX - 1 downto 0);
  signal results_re       : std_logic_vector(C_MAX_RADIX - 1 downto 0);
  signal results_rd_valid : std_logic_vector(C_MAX_RADIX - 1 downto 0);
  signal results_full     : std_logic_vector(C_MAX_RADIX - 1 downto 0);

  -- the preadder always has 5 ports, unused ones are tied off
  type t_preadder_signals_array is array (0 to 4) of t_cmplx;
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

  signal output_mux_out : t_cmplx;   -- selected result-FIFO head

  -- input side counters
  signal phase : std_logic_vector(clogb2(C_MAX_RADIX) - 1 downto 0);
  signal delay_cnt : std_logic_vector(clogb2(C_DELAY_CNT) - 1 downto 0);
  -- output side counters
  signal out_phase : std_logic_vector(clogb2(C_MAX_RADIX) - 1 downto 0);
  signal out_delay_cnt : std_logic_vector(clogb2(C_DELAY_CNT) - 1 downto 0);

  signal twiddle : t_cmplx_twiddle;

  -- flow control: fill phases stall on delay FIFO space, the joint phase on
  -- result FIFO credit. the output side drains X0..X_{r-1} in order.
  signal in_last     : std_logic;
  signal o_ready_int : std_logic;
  signal in_beat     : std_logic;   -- input sample accepted
  signal core_beat   : std_logic;   -- in_beat seen by the core ('0' in bypass)
  signal joint_beat  : std_logic;   -- last-phase beat: butterfly issued
  signal o_valid_int : std_logic;
  signal drain_beat  : std_logic;   -- stored result popped to the output

  -- result FIFO credit: joint beats issued minus last-block pops, counts in-flight words too
  signal res_outstanding : integer range 0 to C_RESULT_DEPTH;
  signal res_credit_ok   : std_logic;

  -- rotator stream: a rot beat launches one word and its twiddle exponent,
  -- the word is delayed C_TW_LAT cycles to meet the twiddle
  signal rot_beat    : std_logic;
  signal rot_word    : t_cmplx;
  signal rot_exp     : integer range 0 to C_MAX_SIZE - 1;
  signal rot_exp_slv : std_logic_vector(C_TW_K_W - 1 downto 0);
  signal rot_out_valid : std_logic;
  signal rot_out_data  : t_cmplx;
  -- twiddle alignment delay line
  type t_rot_align is array (1 to C_TW_LAT) of t_cmplx;
  signal rot_word_d : t_rot_align;
  signal rot_v_d    : std_logic_vector(1 to C_TW_LAT);

  -- what goes into the output skid
  signal skid_we   : std_logic;
  signal skid_data : t_cmplx;

  -- skid credit: launched minus delivered, so the skid can never overflow
  signal out_credit    : integer range 0 to C_SKID_DEPTH;
  signal out_credit_ok : std_logic;

  -- valid pipeline matching the preadder latency
  signal pipe_valid     : std_logic_vector(C_PRE_LAT downto 0);
  signal pipe_valid_out : std_logic;

  signal radix_mask  : std_logic_vector(C_NUM_FIFOS - 1 downto 0);
  signal fifos_we_g  : std_logic_vector(C_NUM_FIFOS - 1 downto 0);
  signal fifos_re_g  : std_logic_vector(C_NUM_FIFOS - 1 downto 0);

  -- radix seen by the output side counters
  signal out_cnt_radix : std_logic_vector(clogb2(C_MAX_RADIX) - 1 downto 0);
begin

  in_last <= '1' when unsigned(phase) = radix_m1 else '0';

  o_ready_int <= o_ready_pre and cfg_ok;

  o_ready <= o_ready_int;
  o_valid <= o_valid_int;

  in_beat    <= i_valid and o_ready_int;
  -- in bypass the core sees no beats, the sample goes straight to the skid
  core_beat  <= in_beat and not bypass;
  joint_beat <= core_beat and in_last;

  -- rotator: word times twiddle from the octant ROM, then a small output skid.
  -- everything on the output is driven from registered state.
  out_credit_ok <= '1' when out_credit < C_SKID_DEPTH else '0';

  PROC_OUT_CREDIT: process(i_clk)
    variable v_inc, v_dec : boolean;
  begin
    if rising_edge(i_clk) then
      if i_reset = '1' then
        out_credit <= 0;
      else
        -- launch: a rot beat, or a direct bypass write
        v_inc := rot_beat = '1' or (bypass = '1' and in_beat = '1');
        v_dec := o_valid_int = '1' and i_ready = '1';
        if v_inc and not v_dec then
          out_credit <= out_credit + 1;
        elsif v_dec and not v_inc then
          out_credit <= out_credit - 1;
        end if;
      end if;
    end if;
  end process PROC_OUT_CREDIT;

  GEN_ROT_MULT: if C_NEED_ROT generate
  begin
    -- delay word and valid to meet the twiddle
    PROC_ROT_REG: process(i_clk)
    begin
      if rising_edge(i_clk) then
        if i_reset = '1' then
          rot_v_d <= (others => '0');
        else
          rot_v_d(1) <= rot_beat;
          for i in 2 to C_TW_LAT loop
            rot_v_d(i) <= rot_v_d(i - 1);
          end loop;
        end if;
        rot_word_d(1) <= rot_word;
        for i in 2 to C_TW_LAT loop
          rot_word_d(i) <= rot_word_d(i - 1);
        end loop;
      end if;
    end process PROC_ROT_REG;

    ROTATOR_INST: entity work.mr_fft_rotator
      generic map (
        G_PIPELINE => G_PIPELINE
      )
      port map (
        i_clk     => i_clk,
        i_reset   => i_reset,
        i_valid   => rot_v_d(C_TW_LAT),
        i_sample  => rot_word_d(C_TW_LAT),
        i_twiddle => twiddle,
        o_valid   => rot_out_valid,
        o_sample  => rot_out_data
      );
    
  else generate
    -- delay 1 stage: W = 1, write straight into the skid
    rot_out_valid <= rot_beat;
    rot_out_data  <= rot_word;
  end generate GEN_ROT_MULT;

  skid_we <= rot_out_valid when bypass = '0' else in_beat;
  skid_data <= rot_out_data when bypass = '0' else input_sample;

  SKID_INST: entity work.mr_fft_fifo
    generic map (
      G_DEPTH => C_SKID_DEPTH
    )
    port map (
      i_clk       => i_clk,
      i_reset     => i_reset,
      i_wr_en     => skid_we,
      i_wr_sample => skid_data,
      i_rd_en     => i_ready,
      o_rd_sample => output_mux_out,
      o_rd_valid  => o_valid_int,
      o_full      => open
    );

  -- capability 1/2 flow: pipelined preadder, result FIFOs, credit based admission
  GEN_FLOW_CAP12: if G_CAPABILITY /= 0 generate
    signal res_head_valid : std_logic;
    signal res_head_word  : t_cmplx;
    -- twiddle exponent p*k for arm p, position k. restarts at every block, never exceeds size.
    signal tw_exp : integer range 0 to C_MAX_SIZE - 1;
  begin

    -- output counters walk all radix result blocks
    out_cnt_radix <= std_logic_vector(to_unsigned(config.radix, clogb2(C_MAX_RADIX)));

    res_credit_ok <= '1' when res_outstanding < C_RESULT_DEPTH else '0';

    -- fill phases wait on delay FIFO space, the joint phase on result credit. i_ready not involved.
    PROC_O_READY: process(bypass, out_credit_ok, in_last, res_credit_ok, fifos_full, phase)
      variable v_idx : integer;
    begin
      v_idx := to_integer(unsigned(phase));
      if bypass = '1' then
        -- bypass: accept while the skid has credit
        o_ready_pre <= out_credit_ok;
      elsif in_last = '1' then
        o_ready_pre <= res_credit_ok;
      elsif v_idx < C_NUM_FIFOS then
        o_ready_pre <= not fifos_full(v_idx);
      else
        -- only during a reconfig transient
        o_ready_pre <= '0';
      end if;
    end process PROC_O_READY;

    -- head of the current result block
    PROC_OUT_HEAD: process(results_rd_valid, results_data_out, out_phase)
      variable v_idx : integer;
    begin
      v_idx := to_integer(unsigned(out_phase));
      if v_idx < C_MAX_RADIX then
        res_head_valid <= results_rd_valid(v_idx);
        res_head_word  <= results_data_out(v_idx);
      else
        res_head_valid <= '0';
        res_head_word  <= (others => (others => '0'));
      end if;
    end process PROC_OUT_HEAD;

    -- pop the result head into the rotator, skid credit guarantees space
    drain_beat <= res_head_valid and out_credit_ok;

    -- bypass beats must not enter the rotator pipeline
    rot_beat <= drain_beat;
    rot_word <= res_head_word;

    GEN_TW_EXP: if C_NEED_ROT generate
      PROC_TW_EXP: process(i_clk)
      begin
        if rising_edge(i_clk) then
          if i_reset = '1' then
            tw_exp <= 0;
          elsif drain_beat = '1' then
            if unsigned(out_delay_cnt) = delay_last_r then
              tw_exp <= 0;   -- block boundary, restart at W^0
            else
              tw_exp <= tw_exp + to_integer(unsigned(out_phase));
            end if;
          end if;
        end if;
      end process PROC_TW_EXP;

      rot_exp <= tw_exp;
    end generate GEN_TW_EXP;

    PROC_RES_OUTSTANDING: process(i_clk)
      variable v_inc, v_dec : boolean;
    begin
      if rising_edge(i_clk) then
        if i_reset = '1' then
          res_outstanding <= 0;
        else
          v_inc := joint_beat = '1';
          v_dec := drain_beat = '1' and
                   unsigned(out_phase) = radix_m1;
          if v_inc and not v_dec then
            res_outstanding <= res_outstanding + 1;
          elsif v_dec and not v_inc then
            res_outstanding <= res_outstanding - 1;
          end if;
        end if;
      end if;
    end process PROC_RES_OUTSTANDING;

    -- valid pipeline next to the preadder pipeline
    pipe_valid(0) <= joint_beat;
    GEN_PIPE_VALID: if C_PRE_LAT > 0 generate
      PROC_PIPE_VALID: process(i_clk)
      begin
        if rising_edge(i_clk) then
          if i_reset = '1' then
            pipe_valid(C_PRE_LAT downto 1) <= (others => '0');
          else
            pipe_valid(C_PRE_LAT downto 1) <= pipe_valid(C_PRE_LAT - 1 downto 0);
          end if;
        end if;
      end process PROC_PIPE_VALID;
    end generate GEN_PIPE_VALID;
    pipe_valid_out <= pipe_valid(C_PRE_LAT);

    -- delay FIFOs: write on fill beats, pop all used arms on joint beats
    GEN_FIFOS_DATA_IN: for i in 0 to C_NUM_FIFOS - 1 generate
      fifos_data_in(i) <= input_demux_out(i);
    end generate GEN_FIFOS_DATA_IN;
    fifos_we_g <= fifos_we when (core_beat = '1' and in_last = '0') else (others => '0');
    fifos_re_g <= radix_mask when joint_beat = '1' else (others => '0');

    -- result FIFOs: write when results land, pop the current head on drain
    results_we(0) <= pipe_valid_out;
    GEN_RESULTS_WE: for j in 1 to C_MAX_RADIX - 1 generate
      results_we(j) <= pipe_valid_out and radix_mask(j - 1);
    end generate GEN_RESULTS_WE;

    PROC_RESULTS_RE: process(drain_beat, out_phase)
      variable v_idx : integer;
    begin
      v_idx := to_integer(unsigned(out_phase));
      results_re <= (others => '0');
      if drain_beat = '1' and v_idx < C_MAX_RADIX then
        results_re(v_idx) <= '1';
      end if;
    end process PROC_RESULTS_RE;

  end generate GEN_FLOW_CAP12;

  -- capability 0 flow: combinational butterfly, X1 written back into the single
  -- delay FIFO, `pending` alternates x and X1 blocks. X0 and drained X1 go
  -- through the output skid.
  GEN_FLOW_RADIX2: if G_CAPABILITY = 0 generate
    signal pending : std_logic;   -- X1 block stored and not yet drained
    -- twiddle exponent for the X1 block, steps of 1
    signal tw_exp  : integer range 0 to C_MAX_SIZE - 1;
  begin

    -- output counters walk the stored X1 block
    out_cnt_radix <= std_logic_vector(radix_m1);
    o_ready_pre <= out_credit_ok when bypass = '1'
                   else (out_credit_ok and not pending) when in_last = '1'
                   else not fifos_full(0);
    drain_beat  <= pending and out_credit_ok;

    -- X0 straight from the butterfly on joint beats, stored X1 on drain beats
    rot_beat <= joint_beat or drain_beat;
    rot_word <= fifos_data_out(0) when pending = '1' else preadder_outputs(0);

    GEN_TW_EXP: if C_NEED_ROT generate
      PROC_TW_EXP: process(i_clk)
      begin
        if rising_edge(i_clk) then
          if i_reset = '1' then
            tw_exp <= 0;
          elsif drain_beat = '1' then
            if unsigned(out_delay_cnt) = delay_last_r then
              tw_exp <= 0;   -- block boundary, restart at W^0
            else
              tw_exp <= tw_exp + 1;
            end if;
          end if;
        end if;
      end process PROC_TW_EXP;

      rot_exp <= tw_exp when pending = '1' else 0;
    end generate GEN_TW_EXP;

    PROC_PENDING: process(i_clk)
    begin
      if rising_edge(i_clk) then
        if i_reset = '1' then
          pending <= '0';
        elsif joint_beat = '1' and
              unsigned(delay_cnt) = delay_last_r then
          pending <= '1';   -- butterfly frame complete, X1 block stored
        elsif drain_beat = '1' and
              unsigned(out_delay_cnt) = delay_last_r then
          pending <= '0';   -- X1 block drained
        end if;
      end if;
    end process PROC_PENDING;

    -- FIFO takes input samples during fill and X1 during the joint phase
    fifos_we_g(0)    <= core_beat;
    fifos_data_in(0) <= preadder_outputs(1) when in_last = '1' else input_demux_out(0);
    fifos_re_g(0)    <= joint_beat or drain_beat;

  end generate GEN_FLOW_RADIX2;

  input_sample <= i_sample;

	GEN_FIFOS: for i in 0 to C_NUM_FIFOS - 1 generate
    FIFO_INST: entity work.mr_fft_fifo
      generic map (
        G_DEPTH      => C_FIFO_DEPTH,
        G_RAM_STYLE  => G_FIFO_RAM_STYLE
      )
      port map (
        i_clk     => i_clk,
        i_reset   => i_reset,
        i_wr_en   => fifos_we_g(i),
        i_wr_sample => fifos_data_in(i),
        i_rd_en   => fifos_re_g(i),
        o_rd_sample => fifos_data_out(i),
        o_rd_valid => fifos_rd_valid(i),
        o_full => fifos_full(i)
      );

  end generate GEN_FIFOS;

  -- result FIFOs exist only in the capability 1/2 flow
  GEN_RESULT_FIFOS_EN: if G_CAPABILITY /= 0 generate
  begin
    GEN_RESULT_FIFOS: for j in 0 to C_MAX_RADIX - 1 generate
      RESULT_FIFO_INST: entity work.mr_fft_fifo
        generic map (
          G_DEPTH      => C_RESULT_DEPTH,
          G_RAM_STYLE  => G_FIFO_RAM_STYLE
        )
        port map (
          i_clk     => i_clk,
          i_reset   => i_reset,
          i_wr_en   => results_we(j),
          i_wr_sample => preadder_outputs(j),
          i_rd_en   => results_re(j),
          o_rd_sample => results_data_out(j),
          o_rd_valid => results_rd_valid(j),
          o_full => results_full(j)
        );

    end generate GEN_RESULT_FIFOS;
  end generate GEN_RESULT_FIFOS_EN;

  -- input demux

  GEN_INPUT_DEMUX_235: if G_CAPABILITY = 2 generate
    PROC_INPUT_DEMUX: process(input_sample, input_demux_sel)
    begin
      input_demux_out <= (others => (others => (others => '0'))); -- Default assignment
      case input_demux_sel is
        when "00" =>
          input_demux_out(0) <= input_sample;
        when "01" =>
          input_demux_out(1) <= input_sample;
        when "10" =>
          input_demux_out(2) <= input_sample;
        when "11" =>
          input_demux_out(3) <= input_sample;
        when others =>
          input_demux_out <= (others => (others => (others => '0'))); -- Default assignment
      end case;
    end process PROC_INPUT_DEMUX;
  end generate GEN_INPUT_DEMUX_235;

  GEN_INPUT_DEMUX_23: if G_CAPABILITY = 1 generate
    PROC_INPUT_DEMUX: process(input_sample, input_demux_sel)
    begin
      input_demux_out <= (others => (others => (others => '0'))); -- Default assignment
      case input_demux_sel is
        when "0" =>
          input_demux_out(0) <= input_sample;
        when "1" =>
          input_demux_out(1) <= input_sample;
        when others =>
          input_demux_out <= (others => (others => (others => '0'))); -- Default assignment
      end case;
    end process PROC_INPUT_DEMUX;
  end generate GEN_INPUT_DEMUX_23;

  GEN_INPUT_DEMUX_2: if G_CAPABILITY = 0 generate
    input_demux_out(0) <= input_sample;
  end generate GEN_INPUT_DEMUX_2;
  
  GEN_PREADDER_INPUTS_235: if G_CAPABILITY = 2 generate
    preadder_s0 <= '1' when s0 = "00" else '0';
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
    preadder_inputs(3) <= (others => (others => '0'));
    preadder_inputs(4) <= (others => (others => '0'));
  end generate GEN_PREADDER_INPUTS_23;

  GEN_PREADDER_INPUTS_2: if G_CAPABILITY = 0 generate
    preadder_inputs(0) <= fifos_data_out(0);
    preadder_inputs(1) <= input_sample;
    preadder_inputs(2) <= (others => (others => '0'));
    preadder_inputs(3) <= (others => (others => '0'));
    preadder_inputs(4) <= (others => (others => '0'));
  end generate GEN_PREADDER_INPUTS_2;

  PREADDER_INST: entity work.mr_fft_preadder
    generic map (
      G_CAPABILITY => G_CAPABILITY,
      G_PIPELINE   => C_PIPELINE
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

  -- registered config decode, per-beat paths start from flops
  PROC_CONFIG_REG: process(i_clk)
  begin
    if rising_edge(i_clk) then
      if i_reset = '1' then
        config_sel <= (others => '0');
        config     <= G_CONFIGS(G_CONFIGS'low);
        bypass     <= '0';
        config_delay <= (others => '0');
        radix_m1   <= to_unsigned(G_CONFIGS(G_CONFIGS'low).radix - 1, radix_m1'length);
      else
        config_sel <= i_config_sel;
        config     <= G_CONFIGS(G_CONFIGS'low + cfg_idx);
        bypass     <= '1' when config.radix < 2 else '0';
        config_delay <= C_DELAY_TBL(cfg_idx);
        radix_m1   <= to_unsigned(config.radix - 1, radix_m1'length);
      end if;
    end if;
  end process PROC_CONFIG_REG;

  -- clamp defends against an out of range sel
  cfg_idx      <= minimum(to_integer(unsigned(config_sel)), C_NUM_CFGS - 1);

  PROC_DELAY_LAST: process(i_clk)
  begin
    if rising_edge(i_clk) then
      delay_last_r <= unsigned(config_delay) - 1;
    end if;
  end process PROC_DELAY_LAST;

  -- count stable cycles after a sel change
  PROC_CFG_SETTLE: process(i_clk)
  begin
    if rising_edge(i_clk) then
      if i_reset = '1' or i_config_sel /= config_sel then
        cfg_settled <= (others => '0');
      else
        cfg_settled <= cfg_settled(C_CFG_SETTLE - 2 downto 0) & '1';
      end if;
    end if;
  end process PROC_CFG_SETTLE;

  -- same-cycle compare, a registered guard alone would be one cycle late
  cfg_ok <= cfg_settled(C_CFG_SETTLE - 1) when i_config_sel = config_sel else '0';

  -- input side counters, advance on accepted samples
  IN_PHASE_DELAY_GEN_INST: entity work.mr_fft_phase_delay_gen
    generic map (
      G_CAPABILITY => G_CAPABILITY,
      G_DELAY_CNT  => C_DELAY_CNT,
      G_MAX_RADIX  => C_MAX_RADIX
    )
    port map (
      i_clk => i_clk,
      i_reset => i_reset,

      i_config_radix => std_logic_vector(to_unsigned(config.radix, clogb2(C_MAX_RADIX))),
      i_config_delay => config_delay,
      i_en => core_beat,

      o_phase => phase,
      o_delay_cnt => delay_cnt
    );

  -- output side counters, advance on drained results
  OUT_PHASE_DELAY_GEN_INST: entity work.mr_fft_phase_delay_gen
    generic map (
      G_CAPABILITY => G_CAPABILITY,
      G_DELAY_CNT  => C_DELAY_CNT,
      G_MAX_RADIX  => C_MAX_RADIX
    )
    port map (
      i_clk => i_clk,
      i_reset => i_reset,

      i_config_radix => out_cnt_radix,
      i_config_delay => config_delay,
      i_en => drain_beat,

      o_phase => out_phase,
      o_delay_cnt => out_delay_cnt
    );

  -- octant twiddle ROM, registered read. not needed in delay 1 stages.
  GEN_TWIDDLE_ROM: if C_NEED_ROT generate
    rot_exp_slv <= std_logic_vector(to_unsigned(rot_exp, C_TW_K_W));

    TWIDDLE_ROM_INST: entity work.twiddle_rom
      generic map (
        CONFIGS => G_CONFIGS,
        SEL_WIDTH => clogb2(C_NUM_CONFIGS),
        K_WIDTH => C_TW_K_W,
        REGISTERED => true,
        G_OUTPUT_REG => C_TW_OUT_REG,
        G_ROM_STYLE => G_TWIDDLE_ROM_STYLE
      )
      port map (
        i_clk => i_clk,

        i_config_sel => config_sel,
        i_k          => rot_exp_slv,

        o_twiddle    => twiddle
      );
  end generate GEN_TWIDDLE_ROM;

  CONTROL_INST: entity work.mr_fft_control
    generic map (
      G_CAPABILITY => G_CAPABILITY,
      G_MAX_RADIX  => C_MAX_RADIX
    )
    port map (
      i_clk => i_clk,
      i_reset => i_reset,

      i_config => config,

      i_phase => phase,

      o_config_s0 => s0,
      o_config_s1 => s1,

      o_input_demux_sel => input_demux_sel,

      o_fifos_we => fifos_we,
      o_radix_mask => radix_mask
    );

  o_sample <= output_mux_out;

  
end architecture rtl;