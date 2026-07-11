library ieee;
use ieee.std_logic_1164.all;
use ieee.fixed_pkg.all;
use ieee.fixed_float_types.all;

library work;
use work.mr_fft_pkg.all;

entity tb_mr_fft_preadder is
end entity tb_mr_fft_preadder;

architecture tb of tb_mr_fft_preadder is
	signal i_x0    : t_cmplx := (re => to_sfixed(0, c_fxp_int_width-1, -c_fxp_frac_width), im => to_sfixed(0, c_fxp_int_width-1, -c_fxp_frac_width));
	signal i_x1    : t_cmplx := (re => to_sfixed(0, c_fxp_int_width-1, -c_fxp_frac_width), im => to_sfixed(0, c_fxp_int_width-1, -c_fxp_frac_width));
	signal i_x2    : t_cmplx := (re => to_sfixed(0, c_fxp_int_width-1, -c_fxp_frac_width), im => to_sfixed(0, c_fxp_int_width-1, -c_fxp_frac_width));
	signal i_x3    : t_cmplx := (re => to_sfixed(0, c_fxp_int_width-1, -c_fxp_frac_width), im => to_sfixed(0, c_fxp_int_width-1, -c_fxp_frac_width));
	signal i_x4    : t_cmplx := (re => to_sfixed(0, c_fxp_int_width-1, -c_fxp_frac_width), im => to_sfixed(0, c_fxp_int_width-1, -c_fxp_frac_width));
	
	signal i_s0    : std_logic_vector(1 downto 0) := "00";
	signal i_s1    : std_logic := '0';

	signal o_X0    : t_cmplx;
	signal o_X1    : t_cmplx;
	signal o_X2    : t_cmplx;
	signal o_X3    : t_cmplx;
	signal o_X4    : t_cmplx;
begin
	dut: entity work.mr_fft_preadder
		generic map (
			G_CAPABILITY => 0
		)
		port map (
			i_x0 => i_x0,
			i_x1 => i_x1,
			i_x2 => i_x2,
			i_x3 => i_x3,
			i_x4 => i_x4,
			i_s0 => i_s0,
			i_s1 => i_s1,
			o_X0 => o_X0,
			o_X1 => o_X1,
			o_X2 => o_X2,
			o_X3 => o_X3,
			o_X4 => o_X4
		);

	stimulus: process
	begin
		i_x0 <= (re => to_sfixed(1.0, c_fxp_int_width-1, -c_fxp_frac_width, fixed_wrap, fixed_truncate), im => to_sfixed(1.0, c_fxp_int_width-1, -c_fxp_frac_width, fixed_wrap, fixed_truncate));
		i_x1 <= (re => to_sfixed(0.5, c_fxp_int_width-1, -c_fxp_frac_width, fixed_wrap, fixed_truncate), im => to_sfixed(0.5, c_fxp_int_width-1, -c_fxp_frac_width, fixed_wrap, fixed_truncate));
		i_x2 <= (re => to_sfixed(0.0, c_fxp_int_width-1, -c_fxp_frac_width, fixed_wrap, fixed_truncate), im => to_sfixed(0.0, c_fxp_int_width-1, -c_fxp_frac_width, fixed_wrap, fixed_truncate));
		i_x3 <= (re => to_sfixed(0.2, c_fxp_int_width-1, -c_fxp_frac_width, fixed_wrap, fixed_truncate), im => to_sfixed(0.2, c_fxp_int_width-1, -c_fxp_frac_width, fixed_wrap, fixed_truncate));
		i_x4 <= (re => to_sfixed(0.1, c_fxp_int_width-1, -c_fxp_frac_width, fixed_wrap, fixed_truncate), im => to_sfixed(0.1, c_fxp_int_width-1, -c_fxp_frac_width, fixed_wrap, fixed_truncate));
		i_s0 <= "00";
		i_s1 <= '1';
		wait for 10 ns;
		wait for 10 ns;

		wait;
	end process stimulus;
end architecture tb;