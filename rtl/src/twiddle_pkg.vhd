library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.fixed_pkg.all;
use ieee.fixed_float_types.all;
use ieee.math_real.all;

library work;
use work.mr_fft_pkg.all;        -- t_cmplx_twiddle, c_twiddle_int_width, c_twiddle_frac_width

-- Twiddle ROM package: builds, at elaboration, one flat ROM with the reduced
-- twiddle set of every config plus per-config base/N/level tables.
-- reduction level: N%4==0 octant (N/8+1 words), even quarter (N/4+1), odd half ((N-1)/2+1).
package twiddle_pkg is

    constant c_twiddle_word_w : natural := 2 * (c_twiddle_int_width + c_twiddle_frac_width);

    type t_coef_rom is array (natural range <>) of t_cmplx_twiddle;
    type t_slv_rom  is array (natural range <>) of std_logic_vector(c_twiddle_word_w-1 downto 0);
    type t_natarr   is array (natural range <>) of natural;

    function f_level(N : natural) return natural;                 -- 3=octant,2=quarter,1=half
    function f_bound(N : natural) return natural;                 -- last stored index
    function f_oct_len(cfgs : t_config_arr) return natural;       -- flat ROM length
    function f_oct_rom(cfgs : t_config_arr) return t_slv_rom;    -- flat coefficient ROM
    function f_oct_base(cfgs : t_config_arr) return t_natarr;     -- per-config base address
    function f_oct_N(cfgs : t_config_arr) return t_natarr;
    function f_oct_level(cfgs : t_config_arr) return t_natarr;
    function f_clog2b(n : natural) return natural;

end package twiddle_pkg;


package body twiddle_pkg is

    constant HI : integer := c_twiddle_int_width - 1;
    constant LO : integer := -c_twiddle_frac_width;

    function f_level(N : natural) return natural is
    begin
        if    N mod 4 = 0 then return 3;
        elsif N mod 2 = 0 then return 2;
        else                   return 1;
        end if;
    end function;

    function f_bound(N : natural) return natural is
    begin
        if    N mod 4 = 0 then return N / 8;
        elsif N mod 2 = 0 then return N / 4;
        else                   return (N - 1) / 2;
        end if;
    end function;

    function f_oct_len(cfgs : t_config_arr) return natural is
        variable s : natural := 0;
    begin
        for i in cfgs'range loop
            s := s + f_bound(cfgs(i).size) + 1;
        end loop;
        return s;
    end function;

    function f_oct_rom(cfgs : t_config_arr) return t_slv_rom is
        variable rom_len: natural := f_oct_len(cfgs);
        variable rom    : t_slv_rom(0 to rom_len-1);
        variable idx    : natural := 0;
        variable N      : natural;
        variable twiddle_re, twiddle_im : std_logic_vector(c_twiddle_word_w/2-1 downto 0);
    begin
        for i in cfgs'range loop
            N := cfgs(i).size;
            for ni in 0 to f_bound(N) loop
                -- round-to-nearest code, then place exactly on the sfixed grid
                twiddle_re := std_logic_vector(to_signed(integer(round( cos(MATH_2_PI * real(ni) / real(N)) * real(2 ** c_twiddle_frac_width))), c_twiddle_word_w/2));
                twiddle_im := std_logic_vector(to_signed(integer(round(-sin(MATH_2_PI * real(ni) / real(N)) * real(2 ** c_twiddle_frac_width))), c_twiddle_word_w/2));
                rom(idx) := twiddle_re & twiddle_im;
                idx := idx + 1;
            end loop;
        end loop;
        return rom;
    end function;

    function f_oct_base(cfgs : t_config_arr) return t_natarr is
        variable b   : t_natarr(cfgs'range);
        variable idx : natural := 0;
    begin
        for i in cfgs'range loop
            b(i) := idx;
            idx  := idx + f_bound(cfgs(i).size) + 1;
        end loop;
        return b;
    end function;

    function f_oct_N(cfgs : t_config_arr) return t_natarr is
        variable a : t_natarr(cfgs'range);
    begin
        for i in cfgs'range loop a(i) := cfgs(i).size; end loop;
        return a;
    end function;

    function f_oct_level(cfgs : t_config_arr) return t_natarr is
        variable a : t_natarr(cfgs'range);
    begin
        for i in cfgs'range loop a(i) := f_level(cfgs(i).size); end loop;
        return a;
    end function;

    function f_clog2b(n : natural) return natural is
        variable r : natural := 0;
        variable v : natural := 1;
    begin
        while v < n loop v := v * 2; r := r + 1; end loop;
        return r;
    end function;

end package body twiddle_pkg;
