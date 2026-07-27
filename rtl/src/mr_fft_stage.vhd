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
    -- pipeline the preadder (see preadder_latency); the stage compensates
    -- with a valid pipeline of the same depth, results land in result FIFOs.
    -- Ignored for G_CAPABILITY = 0: the radix-2 butterfly is one add/sub, so
    -- that stage stays combinational with X1 written back into the delay
    -- FIFO (memory = 1x delay, no result FIFOs).
    G_PIPELINE : boolean := true;
    -- rom_style for the twiddle table: set "block" on slots with big tables
    -- (Vivado leaves inferred ROMs in LUTs by default), "auto" elsewhere
    G_TWIDDLE_ROM_STYLE : string := "auto"
	);
	port (
		i_clk : in  std_logic;
    i_reset : in  std_logic;

    i_config : in  t_config;

    -- index of the current config within G_CONFIGS (selects the twiddle-ROM
    -- block); switch it together with i_config
    i_config_sel : in std_logic_vector(clogb2(G_CONFIGS'length) - 1 downto 0) := (others => '0');

    -- runtime delay (size/radix) for the current config. NOTE: when the value
    -- equals 2**width it wraps to 0; the -1 compares below still work by
    -- modular arithmetic.
    i_config_delay : in std_logic_vector(clogb2(get_delay_cnt(G_CONFIGS)) - 1 downto 0);

    -- Reconfiguration needs NO reset: after a whole number of frames has been
    -- fed and all outputs have been drained (output beats = input beats), the
    -- stage sits in its reset-equivalent state -- both counter pairs wrapped
    -- to (0,0), pending clear, FIFOs empty -- so (i_config, i_config_delay)
    -- may then change to ANY config. Contract: don't offer the next config's
    -- samples before switching the config, and don't switch while the final
    -- drain is still running.

    -- Input and output flows are fully DECOUPLED (separate phase/delay
    -- counters). Delay FIFOs hold only input samples; the (pipelined)
    -- preadder consumes them at the last (radix-1) input phase and its
    -- results X0..X_{r-1} land in per-arm RESULT FIFOs, which the output
    -- side drains block by block whenever the downstream is ready. The input
    -- side never waits for i_ready -- it stalls only on space (delay FIFO
    -- full during fill phases, result-FIFO credit at the joint phase), so
    -- there is no combinational i_ready->o_ready or i_valid->o_valid path.

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
  constant C_DELAY_CNT       : natural := get_delay_cnt(G_CONFIGS);    -- maximum size of the FFT after this stage
  constant C_NUM_CONFIGS     : natural := G_CONFIGS'length; -- Number of configurations
  constant C_FIFO_DEPTH      : natural := C_DELAY_CNT;    -- Depth of each delay FIFO
  -- capability 0 is never pipelined (see G_PIPELINE comment)
  constant C_PIPELINE        : boolean := G_PIPELINE and G_CAPABILITY /= 0;
  -- preadder pipeline depth; the stage's valid pipeline matches it
  constant C_PRE_LAT         : natural := preadder_latency(G_CAPABILITY, C_PIPELINE);
  -- Result FIFOs: one result block (delay words) plus the preadder-pipeline
  -- overlap. Depth D alone would be functionally safe (the credit counter
  -- counts in-flight words, so no overflow), but in back-to-back streaming
  -- the LAST block's final ~C_PRE_LAT+1 pops overlap the next frame's first
  -- landings, so exactly D would stall the joint phase a few cycles per
  -- frame; the +C_PRE_LAT+2 margin makes the schedule stall-free.
  constant C_RESULT_DEPTH    : natural := C_DELAY_CNT + C_PRE_LAT + 2;

  -- A delay-1-only stage (last stage of a chain: C_DELAY_CNT = 1) only ever
  -- rotates by W^0 = 1, so the twiddle ROM, exponent logic and multiplier
  -- are dropped and the rotator degenerates to its data register.
  constant C_NEED_ROT : boolean := C_DELAY_CNT > 1;

  -- Block-RAM twiddle tables get the BRAM output register (fast clock-to-out
  -- instead of the slow latch output); the rotator then aligns with a second
  -- data register. LUT-ROM stages keep the 1-cycle read.
  constant C_TW_OUT_REG : boolean := G_TWIDDLE_ROM_STYLE = "block";
  -- twiddle arrival latency after a rot beat (ROM read + optional output reg)
  -- fold register + address register + BRAM read (+ optional BRAM output reg)
  constant C_TW_LAT     : natural := 3 + boolean'pos(C_TW_OUT_REG);
  -- Rotator pipelining follows the raw G_PIPELINE (NOT C_PIPELINE: capability
  -- 0 keeps its combinational butterfly, but its rotator -- the big radix-2
  -- slots -- is exactly where the DSP pipeline matters).
  constant C_ROT_LAT    : natural := rotator_latency(G_PIPELINE);

  -- rotator: twiddle exponent range and per-radix exponent scaling
  constant C_MAX_SIZE : natural := get_max_size(G_CONFIGS);
  constant C_TW_K_W   : natural := clogb2(C_MAX_SIZE);
  -- output skid after the rotator multiply: its credit gates the beats that
  -- launch words toward the output, counting everything in flight (twiddle
  -- alignment + rotator pipeline), plus margin so streaming never stalls
  constant C_SKID_DEPTH : natural := C_TW_LAT + C_ROT_LAT + 4;

  signal input_sample : t_cmplx;

  signal config_delay : std_logic_vector(clogb2(C_DELAY_CNT) - 1 downto 0);
  -- quasi-static "delay - 1" (block-boundary compare value), registered so
  -- the decrementer sits outside every per-beat comparison
  signal delay_last_r : unsigned(clogb2(C_DELAY_CNT) - 1 downto 0);

  -- delay FIFOs: input samples only (capability 0: X1 write-back as well)
  type t_fifo_data_array is array (0 to C_NUM_FIFOS - 1) of t_cmplx;
  signal fifos_data_in  : t_fifo_data_array;
  signal fifos_data_out : t_fifo_data_array;   -- show-ahead FIFO heads
  signal fifos_we       : std_logic_vector(C_NUM_FIFOS - 1 downto 0);
  signal fifos_rd_valid : std_logic_vector(C_NUM_FIFOS - 1 downto 0);
  signal fifos_full     : std_logic_vector(C_NUM_FIFOS - 1 downto 0);

  -- result FIFOs: butterfly outputs X0..X_{C_MAX_RADIX-1}, one per arm,
  -- written at the preadder pipeline output, drained by the output side
  type t_result_data_array is array (0 to C_MAX_RADIX - 1) of t_cmplx;
  signal results_data_out : t_result_data_array;   -- show-ahead heads
  signal results_we       : std_logic_vector(C_MAX_RADIX - 1 downto 0);
  signal results_re       : std_logic_vector(C_MAX_RADIX - 1 downto 0);
  signal results_rd_valid : std_logic_vector(C_MAX_RADIX - 1 downto 0);
  signal results_full     : std_logic_vector(C_MAX_RADIX - 1 downto 0);

  -- always 5 entries: the preadder entity has 5 input/output ports regardless
  -- of G_CAPABILITY (the unused ones are tied off / left unconnected inside it)
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

  -- INPUT-side counters (drive the demux, FIFO writes and the preadder phase)
  signal phase : std_logic_vector(clogb2(C_MAX_RADIX) - 1 downto 0);
  signal delay_cnt : std_logic_vector(clogb2(C_DELAY_CNT) - 1 downto 0);
  -- OUTPUT-side counters (drive the result drain: phases 0..radix-2)
  signal out_phase : std_logic_vector(clogb2(C_MAX_RADIX) - 1 downto 0);
  signal out_delay_cnt : std_logic_vector(clogb2(C_DELAY_CNT) - 1 downto 0);

  signal twiddle : t_cmplx_twiddle;

  -- Decoupled flow control. Input beats advance the input counters; fill
  -- phases stall only on delay-FIFO space, the joint (last) phase stalls
  -- only on result-FIFO credit. The output side drains result blocks
  -- X0..X_{r-1} in order, throttled purely by head-valid and i_ready.
  signal in_last     : std_logic;
  signal o_ready_int : std_logic;
  signal in_beat     : std_logic;   -- input sample accepted
  signal joint_beat  : std_logic;   -- last-phase beat: butterfly issued
  signal o_valid_int : std_logic;
  signal drain_beat  : std_logic;   -- stored result popped to the output

  -- result-FIFO credit: joint beats issued minus pops of the LAST result
  -- block (an exact upper bound of every result FIFO's occupancy, since the
  -- last block drains last); counts pipeline-in-flight words too.
  signal res_outstanding : integer range 0 to C_RESULT_DEPTH;
  signal res_credit_ok   : std_logic;

  -- rotator stream (driven per flow): a rot beat launches one pre-rotation
  -- word together with its twiddle exponent (address into the octant ROM);
  -- the word is delayed C_TW_LAT cycles to meet the twiddle at the rotator
  signal rot_beat    : std_logic;
  signal rot_word    : t_cmplx;
  signal rot_exp     : integer range 0 to C_MAX_SIZE - 1;
  signal rot_exp_slv : std_logic_vector(C_TW_K_W - 1 downto 0);
  -- twiddle-alignment shift register (depth = the ROM's read latency)
  type t_rot_align is array (1 to C_TW_LAT) of t_cmplx;
  signal rot_word_d : t_rot_align;
  signal rot_v_d    : std_logic_vector(1 to C_TW_LAT);

  -- what actually enters the output skid (rotated word, or the raw word in
  -- a rotator-less delay-1-only stage)
  signal skid_we   : std_logic;
  signal skid_data : t_cmplx;

  -- skid credit: rot beats launched minus words delivered downstream; gating
  -- the beats on it makes skid overflow impossible (in-flight words counted)
  signal out_credit    : integer range 0 to C_SKID_DEPTH;
  signal out_credit_ok : std_logic;

  -- valid pipeline matching the preadder latency: a '1' emerging here means
  -- the preadder outputs carry the results of a joint beat issued C_PRE_LAT
  -- cycles ago
  signal pipe_valid     : std_logic_vector(C_PRE_LAT downto 0);
  signal pipe_valid_out : std_logic;

  signal radix_mask  : std_logic_vector(C_NUM_FIFOS - 1 downto 0);
  signal fifos_we_g  : std_logic_vector(C_NUM_FIFOS - 1 downto 0);
  signal fifos_re_g  : std_logic_vector(C_NUM_FIFOS - 1 downto 0);

  -- radix value fed to the output-side counters (differs per flow scheme)
  signal out_cnt_radix : std_logic_vector(clogb2(C_MAX_RADIX) - 1 downto 0);
begin

  in_last <= '1' when to_integer(unsigned(phase)) = i_config.radix - 1 else '0';

  o_ready <= o_ready_int;
  o_valid <= o_valid_int;

  in_beat    <= i_valid and o_ready_int;
  joint_beat <= in_beat and in_last;

  -- ==================================================================
  -- Rotator (common to both flows): registered pre-rotation word times
  -- the octant-ROM twiddle (1-cycle registered read, addressed by the
  -- rot beat's exponent), landing in a small output skid. o_valid,
  -- o_sample and all gating are functions of registered state.
  -- ==================================================================
  out_credit_ok <= '1' when out_credit < C_SKID_DEPTH else '0';

  PROC_OUT_CREDIT: process(i_clk)
    variable v_inc, v_dec : boolean;
  begin
    if rising_edge(i_clk) then
      if i_reset = '1' then
        out_credit <= 0;
      else
        v_inc := rot_beat = '1';
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
    -- word/valid delayed to meet the twiddle at the rotator (fold register +
    -- ROM read + optional BRAM output register = C_TW_LAT cycles)
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
        o_valid   => skid_we,
        o_sample  => skid_data
      );
  else generate
    -- delay-1-only stage: W = 1 -- no ROM, no rotator, no data register;
    -- the beat writes the word directly into the skid
    skid_we   <= rot_beat;
    skid_data <= rot_word;
  end generate GEN_ROT_MULT;

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

  -- ==================================================================
  -- Capability 1/2 flow (radix23, radix235): pipelined preadder, results land
  -- in per-arm result FIFOs, credit-based admission, head-valid drain
  -- ==================================================================
  GEN_FLOW_CAP12: if G_CAPABILITY /= 0 generate
    signal res_head_valid : std_logic;
    signal res_head_word  : t_cmplx;
    -- Twiddle exponent p*k for output beat (arm p, position k): steps by the
    -- arm index within a block, restarts at 0 at each block start. The config
    -- tables carry per-slot sub-sizes, so p*k < size always -- no modulo.
    signal tw_exp : integer range 0 to C_MAX_SIZE - 1;
  begin

    -- output-side counters walk all radix result blocks X0..X_{r-1}
    out_cnt_radix <= std_logic_vector(to_unsigned(i_config.radix, clogb2(C_MAX_RADIX)));

    res_credit_ok <= '1' when res_outstanding < C_RESULT_DEPTH else '0';

    -- fill phases: accept while the target delay FIFO has space; joint
    -- phase: accept while the result FIFOs have credit. i_ready is NOT
    -- involved.
    PROC_O_READY: process(in_last, res_credit_ok, fifos_full, phase)
      variable v_idx : integer;
    begin
      v_idx := to_integer(unsigned(phase));
      if in_last = '1' then
        o_ready_int <= res_credit_ok;
      elsif v_idx < C_NUM_FIFOS then
        o_ready_int <= not fifos_full(v_idx);
      else
        -- reconfig transient only: new (smaller) radix applied while the old
        -- phase value is still in the counter
        o_ready_int <= '0';
      end if;
    end process PROC_O_READY;

    -- output side: present the head of the current result block
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

    -- a drain beat pops the result head into the rotator (skid credit
    -- guarantees it has somewhere to land; i_ready only drains the skid)
    drain_beat <= res_head_valid and out_credit_ok;

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
              tw_exp <= 0;   -- block boundary: next block restarts at W^0
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
                   to_integer(unsigned(out_phase)) = i_config.radix - 1;
          if v_inc and not v_dec then
            res_outstanding <= res_outstanding + 1;
          elsif v_dec and not v_inc then
            res_outstanding <= res_outstanding - 1;
          end if;
        end if;
      end if;
    end process PROC_RES_OUTSTANDING;

    -- valid pipeline alongside the (free-running) preadder pipeline
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

    -- delay-FIFO side: input samples only; writes on fill beats, pops (all
    -- used arms) on joint beats; the drain never touches the delay FIFOs
    GEN_FIFOS_DATA_IN: for i in 0 to C_NUM_FIFOS - 1 generate
      fifos_data_in(i) <= input_demux_out(i);
    end generate GEN_FIFOS_DATA_IN;
    fifos_we_g <= fifos_we when (in_beat = '1' and in_last = '0') else (others => '0');
    fifos_re_g <= radix_mask when joint_beat = '1' else (others => '0');

    -- result-FIFO enables: writes when the preadder results land (X0 always,
    -- X1..X_{r-1} per the radix), pop of the current block's head on drain
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

  -- ==================================================================
  -- Radix-2-only flow (capability 0): as before pipelining, plus an
  -- output skid. The one-add butterfly stays combinational; X1 is
  -- written back into the single delay FIFO (memory = 1x delay, no
  -- result FIFOs) and the `pending` flop interlocks the frames (x and
  -- X1 blocks alternate in the FIFO). X0 (joint beats) and drained X1
  -- words feed a 2-deep show-ahead skid FIFO whose registered state
  -- drives o_valid/o_sample and gates the joint phase -- so there is no
  -- combinational i_valid->o_valid or i_ready->o_ready path here either.
  -- ==================================================================
  GEN_FLOW_RADIX2: if G_CAPABILITY = 0 generate
    signal pending : std_logic;   -- X1 block stored and not yet drained
    -- twiddle exponent for the X1 block (arm 1): k, steps of 1, k < size
    signal tw_exp  : integer range 0 to C_MAX_SIZE - 1;
  begin

    -- output-side counters walk the single stored block (X1)
    out_cnt_radix <= std_logic_vector(to_unsigned(i_config.radix - 1, clogb2(C_MAX_RADIX)));

    o_ready_int <= (out_credit_ok and not pending) when in_last = '1'
                   else not fifos_full(0);
    drain_beat  <= pending and out_credit_ok;

    -- rotator stream: X0 straight from the butterfly on joint beats
    -- (exponent 0 -> W = 1), stored X1 words on drain beats (exclusive
    -- via pending)
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
              tw_exp <= 0;   -- block boundary: next block restarts at W^0
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

    -- the FIFO takes the input sample during the fill phase and the X1
    -- write-back during the joint phase; joint pops (preadder input) and
    -- drain pops (into the skid) are mutually exclusive via pending
    fifos_we_g(0)    <= in_beat;
    fifos_data_in(0) <= preadder_outputs(1) when in_last = '1' else input_demux_out(0);
    fifos_re_g(0)    <= joint_beat or drain_beat;

  end generate GEN_FLOW_RADIX2;

  input_sample <= i_sample;

	GEN_FIFOS: for i in 0 to C_NUM_FIFOS - 1 generate
    FIFO_INST: entity work.mr_fft_fifo
      generic map (
        G_DEPTH      => C_FIFO_DEPTH
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
          G_DEPTH      => C_RESULT_DEPTH
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

  -- Input demux logic

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

  config_delay <= i_config_delay;

  PROC_DELAY_LAST: process(i_clk)
  begin
    if rising_edge(i_clk) then
      delay_last_r <= unsigned(config_delay) - 1;
    end if;
  end process PROC_DELAY_LAST;

  -- input-side counters: advance on accepted input samples
  IN_PHASE_DELAY_GEN_INST: entity work.mr_fft_phase_delay_gen
    generic map (
      G_CAPABILITY => G_CAPABILITY,
      G_DELAY_CNT  => C_DELAY_CNT,
      G_MAX_RADIX  => C_MAX_RADIX
    )
    port map (
      i_clk => i_clk,
      i_reset => i_reset,

      i_config_radix => std_logic_vector(to_unsigned(i_config.radix, clogb2(C_MAX_RADIX))),
      i_config_delay => config_delay,
      i_en => in_beat,

      o_phase => phase,
      o_delay_cnt => delay_cnt
    );

  -- output-side counters: advance on drained results (capability 1/2: all radix
  -- result blocks X0..X_{r-1}; radix-2 flow: the single stored X1 block)
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

  -- octant twiddle ROM: addressed by the rot beat's exponent; REGISTERED
  -- read aligns with the rotator's data register (multiply happens one
  -- cycle after the beat). Absent in delay-1-only stages.
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

        i_config_sel => i_config_sel,
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

      i_config => i_config,

      i_phase => phase,

      o_config_s0 => s0,
      o_config_s1 => s1,

      o_input_demux_sel => input_demux_sel,

      o_fifos_we => fifos_we,
      o_radix_mask => radix_mask
    );

  o_sample <= output_mux_out;

  
end architecture rtl;