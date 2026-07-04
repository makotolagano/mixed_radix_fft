library ieee;
use ieee.std_logic_1164.all;

library work;
use work.mr_fft_pkg.all;

entity mr_fft_preadder is
    port (
        i_x0    : in  t_cmplx;
        i_x1    : in  t_cmplx;
        i_x2    : in  t_cmplx;
        i_x3    : in  t_cmplx;
        i_x4    : in  t_cmplx;

        i_s0    : in  std_logic_vector(1 downto 0);
        i_s1    : in  std_logic;

        o_X0    : out t_cmplx;
        o_X1    : out t_cmplx;
        o_X2    : out t_cmplx;
        o_X3    : out t_cmplx;
        o_X4    : out t_cmplx
    );
end entity mr_fft_preadder;

architecture rtl of mr_fft_preadder is
begin
    -- Some example t_cmplx arithemtic to check if fixed point complex operations work

    o_X0 <= i_x0 + i_x1;
    o_X1 <= i_x0 - i_x1;

end architecture rtl;