library ieee;
use ieee.std_logic_1164.all;
use ieee.fixed_pkg.all;

-- Custom packages
library work;
use work.mr_fft_pkg.all;


entity mr_fft_rotator is
	port (
		i_sample 	: in  t_cmplx;
		i_twiddle : in  t_cmplx_twiddle;
		o_sample 	: out t_cmplx
	);
end entity mr_fft_rotator;

architecture rtl of mr_fft_rotator is

begin

	o_sample <= i_sample * i_twiddle;

end architecture rtl;