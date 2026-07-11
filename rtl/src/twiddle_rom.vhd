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
-- Pipeline (single register so block RAM infers at 1-cycle latency):
--   stage 0 (comb)          : fold k -> ROM address + carry conj/unit flags
--   stage 1 (read)          : ROM read -- the ONLY register
--                             REGISTERED=true  -> clocked read  -> block RAM, 1-cycle latency
--                             REGISTERED=false -> comb read      -> distributed ROM, 0 latency
--   stage 2 (comb)          : reconstruct from the ROM word; drives the output directly
-- Read latency is at most one clock; the conj/unit flags travel with the ROM
-- word so reconstruct stays aligned in both modes.
-- ---------------------------------------------------------------------------
entity twiddle_rom is
    generic (
        CONFIGS    : t_config_arr := c_stage4_cfgs;               -- all (radix,size) configs for this stage
        SEL_WIDTH  : natural := 6;                    -- config_sel width (>= clog2(CONFIGS'length))
        K_WIDTH    : natural := 12;                    -- k width (>= clog2(max config size))
        REGISTERED : boolean := true
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
    constant C_ROM  : t_slv_rom := f_oct_rom(CONFIGS);
    constant C_BASE : t_natarr  := f_oct_base(CONFIGS);
    constant C_NN   : t_natarr  := f_oct_N(CONFIGS);
    constant C_LVL  : t_natarr  := f_oct_level(CONFIGS);
    constant C_WW   : natural    := i_k'length + 3;          -- fold working width (headroom 8*r)
    constant C_HI   : integer    := c_twiddle_int_width - 1;
    constant C_LO   : integer    := -c_twiddle_frac_width;
    constant C_CW   : natural    := c_twiddle_int_width + c_twiddle_frac_width;  -- bits per component

    -- Block RAM inference
    -- attribute rom_style : string;
    -- attribute rom_style of ROM : constant is "block";

    -- stage 0 -> stage 1 : folded ROM address + reconstruct flags
    signal addr      : natural range 0 to C_ROM'length-1;
    signal conjugate : std_logic;
    signal m         : unsigned(1 downto 0); -- unit j^m

    -- stage 1 (ROM read) outputs -- the single register when REGISTERED = true
    signal rom_q       : std_logic_vector(c_twiddle_word_w-1 downto 0);
    signal conjugate_q : std_logic;
    signal m_q         : unsigned(1 downto 0);

    -- stage 2 (reconstruct) output
    signal twiddle_c   : t_cmplx_twiddle;
begin

    -- ---- stage 0 (comb): fold k -> ROM address, carry conj/unit flags ----
    fold : process (i_config_sel, i_k)
        variable v_sel, v_base_a, v_lv : natural;
        variable v_N, v_N2, v_N4, v_r    : unsigned(C_WW-1 downto 0);
        variable v_m               : unsigned(1 downto 0);         -- unit j^m; +k wraps mod 4
        variable v_conjugate       : std_logic;
    begin
        v_sel    := to_integer(unsigned(i_config_sel));
        v_base_a := C_BASE(v_sel);
        v_lv     := C_LVL(v_sel);
        v_N      := to_unsigned(C_NN(v_sel), C_WW);
        v_N2     := shift_right(v_N, 1);                             -- N/2
        v_N4     := shift_right(v_N, 2);                             -- N/4

        v_r  := resize(unsigned(i_k), C_WW);                           -- k < N by construction
        v_m  := "00";
        v_conjugate := '0';

        if shift_left(v_r, 1) > v_N then               -- conjugate fold
            v_r := v_N - v_r;   v_conjugate := not v_conjugate;
        end if;
        if v_lv >= 2 and shift_left(v_r, 2) > v_N then -- quarter fold (t=-1)
            v_r := v_N2 - v_r;  v_m := v_m + 2;  v_conjugate := not v_conjugate;
        end if;
        if v_lv >= 3 and shift_left(v_r, 3) > v_N then -- octant fold (t=-j; +j if conjugate)
            v_r := v_N4 - v_r;
            if v_conjugate = '1' then v_m := v_m + 1; else v_m := v_m + 3; end if;
            v_conjugate := not v_conjugate;
        end if;

        addr      <= v_base_a + to_integer(v_r);
        conjugate <= v_conjugate;
        m         <= v_m;
    end process;

    -- ---- stage 1: ROM read (the ONLY register).  Registered -> block RAM, 1-cycle latency;
    --      combinational -> distributed ROM, 0 latency.  Flags travel with the ROM word.
    read_reg_g : if REGISTERED generate
        process (i_clk) begin
            if rising_edge(i_clk) then
                rom_q <= C_ROM(addr);
                conjugate_q  <= conjugate;
                m_q   <= m;
            end if;
        end process;
    end generate;

    read_comb_g : if not REGISTERED generate
        rom_q <= C_ROM(addr);
        conjugate_q  <= conjugate;
        m_q   <= m;
    end generate;

    -- ---- stage 2 (comb): unpack the ROM word and reconstruct; drives the output directly ----
    recon : process (rom_q, conjugate_q, m_q)
        variable a, b, na, nb : sfixed(C_HI downto C_LO);
    begin
        a := to_sfixed(rom_q(2*C_CW-1 downto C_CW), C_HI, C_LO); -- re
        b := to_sfixed(rom_q(C_CW-1   downto  0), C_HI, C_LO);   -- im
        if conjugate_q = '1' then -- conjugate: negate Im
            b := resize(-b, C_HI, C_LO, fixed_saturate, fixed_truncate);
        end if;
        na := resize(-a, C_HI, C_LO, fixed_saturate, fixed_truncate);
        nb := resize(-b, C_HI, C_LO, fixed_saturate, fixed_truncate);

        case m_q is                                              -- apply unit j^m: swapping/negation
            when "00"   => twiddle_c.re <= a;  twiddle_c.im <= b;   --  1 : ( a,  b)
            when "01"   => twiddle_c.re <= nb; twiddle_c.im <= a;   --  j : (-b,  a)
            when "10"   => twiddle_c.re <= na; twiddle_c.im <= nb;  -- -1 : (-a, -b)
            when others => twiddle_c.re <= b;  twiddle_c.im <= na;  -- -j : ( b, -a)
        end case;
    end process;

    o_twiddle <= twiddle_c;

end architecture;
