library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.fixed_pkg.all;
use ieee.fixed_float_types.all;

library work;
use work.mr_fft_pkg.all;        -- t_cmplx_twiddle
use work.twiddle_pkg.all;   -- t_config_arr, t_coef_rom, f_oct_*

-- ---------------------------------------------------------------------------
-- Octant-symmetry stored twiddle ROM.
--
-- Holds ALL of a stage's configs (flat t_cmplx_twiddle ROM + per-config base/N/level
-- tables, built at elaboration). Runtime datapath is explicit logic: read per-config
-- params, fold the exponent, read the ROM, reconstruct with conjugate/negate/swap.
--
--   inputs : config_sel (which config), k = phase*delay_cnt*R^stage_index  (0..N-1)
--   fold   : reduce k -> stored index r, flags conj + unit-code m (0..3: unit j^m)
--   read   : ROM[ base(config_sel) + r ]  -- SYNCHRONOUS (registered) so it maps to BRAM
--   recon  : if conj negate Im ; then apply j^m (swap/negate) -- muxes, no multiplier
--   k = 0 (arm 0) -> address base -> W = 1.
--
-- Pipeline (REGISTERED=true; REGISTERED=false leaves every stage comb),
-- one fold per stage so each cycle holds one compare + subtract:
--   stage 0a (comb)  : conjugate fold
--   stage 0b (reg+comb): fold register, then quarter fold
--   stage 0c (reg+comb): fold register, then octant fold + base add
--   stage 0d (reg)   : address register (isolates BRAM address setup)
--   stage 1  (reg)   : ROM read -- clocked -> block RAM
--   stage 1b (reg)   : optional BRAM output register (G_OUTPUT_REG)
--   stage 2  (reg)   : reconstruct (negate/swap muxes), registered output
-- REGISTERED latency: 5 cycles (+1 with G_OUTPUT_REG); the conj/unit flags
-- travel with the ROM word so reconstruct stays aligned in both modes.
-- ---------------------------------------------------------------------------
entity twiddle_rom is
    generic (
        CONFIGS    : t_config_arr;                                -- all (radix,size) configs for this stage
        SEL_WIDTH  : natural := 6;                    -- config_sel width (>= clog2(CONFIGS'length))
        K_WIDTH    : natural := 12;                    -- k width (>= clog2(max config size))
        REGISTERED : boolean := true;
        -- second read register stage: Vivado absorbs it as the BRAM's optional
        -- output register (DOA_REG=1, fast clock-to-out instead of the slow
        -- latch output). +1 cycle read latency; only meaningful with
        -- REGISTERED and a block-RAM table.
        G_OUTPUT_REG : boolean := false;
        -- Vivado rom_style for the coefficient table: inferred ROMs default to
        -- LUTs even when large, so set "block" per instance where the table is
        -- big (radix2/radix23 slots); keep "auto"/"distributed" for small ones
        G_ROM_STYLE : string := "auto"
    );
    port (
        i_clk        : in  std_logic := '0';
        i_config_sel : in  std_logic_vector(SEL_WIDTH-1 downto 0);  -- index into CONFIGS
        i_k          : in  std_logic_vector(K_WIDTH-1 downto 0);    -- twiddle exponent (0 .. N-1)
        o_twiddle    : out t_cmplx_twiddle          -- reconstructed twiddle (sfixed re/im)
    );
end entity;

architecture rtl of twiddle_rom is
    -- Flat coefficient ROM as packed std_logic_vector words (re & im) so a registered
    -- read maps to block RAM. (An array of records does not infer BRAM in Vivado.)
    -- Held in a SIGNAL with an initial value (never written): Vivado does not
    -- reliably map reads of a CONSTANT array onto RAMB primitives. The style
    -- comes from the ROM_STYLE generic (per instance).
    constant C_ROM_INIT : t_slv_rom := f_oct_rom(CONFIGS);
    signal   rom_mem    : t_slv_rom(C_ROM_INIT'range) := C_ROM_INIT;

attribute rom_style : string;
    attribute rom_style of rom_mem : signal is G_ROM_STYLE;
    constant C_BASE : t_natarr  := f_oct_base(CONFIGS);
    constant C_NN   : t_natarr  := f_oct_N(CONFIGS);
    constant C_LVL  : t_natarr  := f_oct_level(CONFIGS);
    constant C_WW   : natural    := i_k'length + 3;          -- fold working width (headroom 8*r)
    constant C_AW   : natural    := maximum(1, f_clog2b(C_ROM_INIT'length));  -- addr bits
    constant C_PW   : natural    := maximum(C_WW, C_AW) + 1; -- addr arithmetic width
    constant C_HI   : integer    := c_twiddle_int_width - 1;
    constant C_LO   : integer    := -c_twiddle_frac_width;
    constant C_CW   : natural    := c_twiddle_int_width + c_twiddle_frac_width;  -- bits per component


    -- per-config fold parameters, registered when REGISTERED (i_config_sel
    -- is quasi-static, so the 53-entry table reads leave the per-beat path)
    signal base_p          : natural range 0 to C_ROM_INIT'length-1;
    signal lvl_p           : natural range 1 to 3;
    signal n_p, n2_p, n4_p : unsigned(C_WW-1 downto 0);
    signal n8_p            : unsigned(C_WW-1 downto 0);
    signal bn4_p           : unsigned(C_PW-1 downto 0);      -- base + N/4 (octant branch)

    -- stage 0a (conjugate fold) -> stage 0b (the conjugate fold never
    -- touches m, so no m here)
    signal r1, r1_s    : unsigned(C_WW-1 downto 0);
    signal conj1, conj1_s : std_logic;

    -- stage 0b (quarter fold) -> stage 0c
    signal r2, r2_s    : unsigned(C_WW-1 downto 0);
    signal conj2, conj2_s : std_logic;
    signal m2, m2_s    : unsigned(1 downto 0);

    -- stage 0c -> stage 0d : folded ROM address + reconstruct flags. The
    -- address wraps (unsigned truncation) during a config-switch transient;
    -- the guarded read below keeps the discarded transient legal in sim.
    signal addr      : unsigned(C_AW-1 downto 0);
    signal conjugate : std_logic;
    signal m         : unsigned(1 downto 0); -- unit j^m

    -- stage 0d (address register) -> stage 1: the read cycle then holds only
    -- the register -> BRAM routing plus the read guard
    signal addr_r : unsigned(C_AW-1 downto 0);
    signal conj_b : std_logic;
    signal m_b    : unsigned(1 downto 0);

    -- stage 1 (ROM read) outputs -- the single register when REGISTERED = true
    signal rom_q       : std_logic_vector(c_twiddle_word_w-1 downto 0);
    signal conjugate_q : std_logic;
    signal m_q         : unsigned(1 downto 0);

    -- optional stage 1b (BRAM output register) -- pass-through when disabled
    signal rom_q2       : std_logic_vector(c_twiddle_word_w-1 downto 0);
    signal conjugate_q2 : std_logic;
    signal m_q2         : unsigned(1 downto 0);

    -- stage 2 (reconstruct) output
    signal twiddle_c   : t_cmplx_twiddle;
begin

    -- ---- per-config fold parameters (table reads off the per-beat path) ----
    params_reg_g : if REGISTERED generate
        process (i_clk)
            variable v_sel : natural;
            variable v_N   : unsigned(C_WW-1 downto 0);
        begin
            if rising_edge(i_clk) then
                v_sel  := to_integer(unsigned(i_config_sel));
                v_N    := to_unsigned(C_NN(v_sel), C_WW);
                base_p <= C_BASE(v_sel);
                lvl_p  <= C_LVL(v_sel);
                n_p    <= v_N;
                n2_p   <= shift_right(v_N, 1);                       -- N/2
                n4_p   <= shift_right(v_N, 2);                       -- N/4
                n8_p   <= shift_right(v_N, 3);                       -- N/8
                bn4_p  <= to_unsigned(C_BASE(v_sel) + C_NN(v_sel) / 4, C_PW);
            end if;
        end process;
    end generate;

    params_comb_g : if not REGISTERED generate
        process (i_config_sel)
            variable v_sel : natural;
            variable v_N   : unsigned(C_WW-1 downto 0);
        begin
            v_sel  := to_integer(unsigned(i_config_sel));
            v_N    := to_unsigned(C_NN(v_sel), C_WW);
            base_p <= C_BASE(v_sel);
            lvl_p  <= C_LVL(v_sel);
            n_p    <= v_N;
            n2_p   <= shift_right(v_N, 1);
            n4_p   <= shift_right(v_N, 2);
            n8_p   <= shift_right(v_N, 3);
            bn4_p  <= to_unsigned(C_BASE(v_sel) + C_NN(v_sel) / 4, C_PW);
        end process;
    end generate;

    -- ---- stage 0a (comb): conjugate fold (2k > N <=> k > N/2 exactly) ----
    fold_a : process (i_k, n_p, n2_p)
        variable v_r : unsigned(C_WW-1 downto 0);
    begin
        v_r := resize(unsigned(i_k), C_WW);                          -- k < N by construction
        if v_r > n2_p then
            r1    <= n_p - v_r;
            conj1 <= '1';
        else
            r1    <= v_r;
            conj1 <= '0';
        end if;
    end process;

    fold_a_reg_g : if REGISTERED generate
        process (i_clk) begin
            if rising_edge(i_clk) then
                r1_s    <= r1;
                conj1_s <= conj1;
            end if;
        end process;
    end generate;

    fold_a_comb_g : if not REGISTERED generate
        r1_s    <= r1;
        conj1_s <= conj1;
    end generate;

    -- ---- stage 0b (comb): quarter fold, t=-1 (4r > N <=> r > N/4) ----
    fold_b : process (r1_s, conj1_s, n2_p, n4_p, lvl_p)
    begin
        if lvl_p >= 2 and r1_s > n4_p then
            r2    <= n2_p - r1_s;
            conj2 <= not conj1_s;
            m2    <= "10";
        else
            r2    <= r1_s;
            conj2 <= conj1_s;
            m2    <= "00";
        end if;
    end process;

    fold_b_reg_g : if REGISTERED generate
        process (i_clk) begin
            if rising_edge(i_clk) then
                r2_s    <= r2;
                conj2_s <= conj2;
                m2_s    <= m2;
            end if;
        end process;
    end generate;

    fold_b_comb_g : if not REGISTERED generate
        r2_s    <= r2;
        conj2_s <= conj2;
        m2_s    <= m2;
    end generate;

    -- ---- stage 0c (comb): octant fold + base add -> ROM address ----
    -- folded address is (base + N/4) - r with bn4_p precomputed, so both
    -- branches are ONE adder deep and run in parallel; the compare
    -- (8r > N <=> r > N/8) only drives the selects
    fold_c : process (r2_s, conj2_s, m2_s, n8_p, lvl_p, base_p, bn4_p)
        variable v_fold : boolean;
    begin
        v_fold := lvl_p >= 3 and r2_s > n8_p;            -- octant fold (t=-j; +j if conjugate)

        if v_fold then
            addr      <= resize(bn4_p - r2_s, C_AW);
            conjugate <= not conj2_s;
            if conj2_s = '1' then m <= m2_s + 1; else m <= m2_s + 3; end if;
        else
            addr      <= resize(to_unsigned(base_p, C_PW) + r2_s, C_AW);
            conjugate <= conj2_s;
            m         <= m2_s;
        end if;
    end process;

    -- ---- stage 0d: address register (REGISTERED mode) -- isolates the BRAM
    --      address setup (or LUT-ROM decode) from the fold logic ----
    addr_reg_g : if REGISTERED generate
        process (i_clk) begin
            if rising_edge(i_clk) then
                addr_r <= addr;
                conj_b <= conjugate;
                m_b    <= m;
            end if;
        end process;
    end generate;

    addr_comb_g : if not REGISTERED generate
        addr_r <= addr;
        conj_b <= conjugate;
        m_b    <= m;
    end generate;

    -- ---- stage 1: ROM read.  Registered -> block RAM, 1-cycle latency;
    --      combinational -> distributed ROM, 0 latency.  Flags travel with the ROM word.
    read_reg_g : if REGISTERED generate
        process (i_clk) begin
            if rising_edge(i_clk) then
                -- read guard: a wrapped transient address (config switch)
                -- reads nothing -- keeps the array index legal in simulation
                if to_integer(addr_r) <= C_ROM_INIT'length - 1 then
                    rom_q <= rom_mem(to_integer(addr_r));
                end if;
                conjugate_q  <= conj_b;
                m_q   <= m_b;
            end if;
        end process;
    end generate;

    read_comb_g : if not REGISTERED generate
        rom_q <= rom_mem(minimum(to_integer(addr_r), C_ROM_INIT'length - 1));
        conjugate_q  <= conj_b;
        m_q   <= m_b;
    end generate;

    -- ---- stage 1b: optional BRAM output register (flags travel along) ----
    out_reg_g : if G_OUTPUT_REG generate
        process (i_clk) begin
            if rising_edge(i_clk) then
                rom_q2       <= rom_q;
                conjugate_q2 <= conjugate_q;
                m_q2         <= m_q;
            end if;
        end process;
    end generate;

    out_comb_g : if not G_OUTPUT_REG generate
        rom_q2       <= rom_q;
        conjugate_q2 <= conjugate_q;
        m_q2         <= m_q;
    end generate;

    -- ---- stage 2: unpack the ROM word and reconstruct ----
    -- negations in parallel off the raw word, then a 4:1 mux per component:
    -- one carry-chain level (chaining conjugate then j^m would need two)
    recon : process (rom_q2, conjugate_q2, m_q2)
        variable a, b, na, nb : sfixed(C_HI downto C_LO);
    begin
        a  := to_sfixed(rom_q2(2*C_CW-1 downto C_CW), C_HI, C_LO); -- re
        b  := to_sfixed(rom_q2(C_CW-1   downto  0), C_HI, C_LO);   -- im
        na := resize(-a, C_HI, C_LO, fixed_saturate, fixed_truncate);
        nb := resize(-b, C_HI, C_LO, fixed_saturate, fixed_truncate);

        -- twiddle = j^m * (a, conj ? -b : b), expanded per component
        case m_q2 is
            when "00" =>                                        --  1 : ( x,  y)
                twiddle_c.re <= a;
                twiddle_c.im <= nb when conjugate_q2 = '1' else b;
            when "01" =>                                        --  j : (-y,  x)
                twiddle_c.re <= b  when conjugate_q2 = '1' else nb;
                twiddle_c.im <= a;
            when "10" =>                                        -- -1 : (-x, -y)
                twiddle_c.re <= na;
                twiddle_c.im <= b  when conjugate_q2 = '1' else nb;
            when others =>                                      -- -j : ( y, -x)
                twiddle_c.re <= nb when conjugate_q2 = '1' else b;
                twiddle_c.im <= na;
        end case;
    end process;

    -- reconstruct output register (REGISTERED mode)
    recon_reg_g : if REGISTERED generate
        process (i_clk) begin
            if rising_edge(i_clk) then
                o_twiddle <= twiddle_c;
            end if;
        end process;
    end generate;

    recon_comb_g : if not REGISTERED generate
        o_twiddle <= twiddle_c;
    end generate;

end architecture;
