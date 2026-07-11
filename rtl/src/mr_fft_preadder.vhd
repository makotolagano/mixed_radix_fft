library ieee;
use ieee.std_logic_1164.all;
use ieee.fixed_pkg.all;

-- Custom packages
library work;
use work.mr_fft_pkg.all;


entity mr_fft_preadder is
	generic (
		G_CAPABILITY : natural := 0
	);
	port (
		i_x0 : in  t_cmplx;
		i_x1 : in  t_cmplx;
		i_x2 : in  t_cmplx;
		i_x3 : in  t_cmplx;
		i_x4 : in  t_cmplx;

		i_s0 : in  std_logic_vector(1 downto 0);
		i_s1 : in  std_logic;

		o_X0 : out t_cmplx;
		o_X1 : out t_cmplx;
		o_X2 : out t_cmplx;
		o_X3 : out t_cmplx;
		o_X4 : out t_cmplx
	);
end entity mr_fft_preadder;

architecture rtl of mr_fft_preadder is
	-- Type definitions
	type t_cmplx_wide_array is array (5 downto 0) of t_cmplx_wide;
	
	signal input_x0_wide, input_x1_wide, input_x2_wide, input_x3_wide, input_x4_wide 			: t_cmplx_wide;
	signal output_X0_wide, output_X1_wide, output_X2_wide, output_X3_wide, output_X4_wide : t_cmplx_wide;

begin

	-- Resize input signals to wide fixed-point representation
	input_x0_wide <= resize(i_x0, input_x0_wide);
	input_x1_wide <= resize(i_x1, input_x1_wide);
	input_x2_wide <= resize(i_x2, input_x2_wide);
	input_x3_wide <= resize(i_x3, input_x3_wide);
	input_x4_wide <= resize(i_x4, input_x4_wide);

	GEN_235: if G_CAPABILITY = 2 generate
		PROC_CALC_235: process(input_x0_wide, input_x1_wide, input_x2_wide, input_x3_wide, input_x4_wide, i_s0, i_s1)
			variable t2_s0_mux_out, t2_mul1_out, t2_mul2_out, t2_mul3_out, t2_mul4_out : t_cmplx_wide;
			variable t0, t1, t2, t3 																									 : t_cmplx_wide_array;
		begin
			-- Stage 0
			t0(0) := input_x0_wide;
			t0(1) := input_x1_wide + input_x4_wide;
			t0(2) := input_x2_wide + input_x3_wide;
			t0(3) := input_x1_wide - input_x4_wide;
			t0(4) := input_x2_wide - input_x3_wide;

			-- Stage 1
			t1(0) := t0(0);
			t1(1) := t0(1) + t0(2);
			t1(2) := t0(1) - t0(2);
			t1(3) := t0(3);
			t1(4) := t0(4);
			t1(5) := t0(3) + t0(4);

			-- Stage 2
			t2(0) := t1(0) + t1(1);
			
			case i_s0 is
				when "00" => -- no shift
					t2_s0_mux_out := t1(1);
				when "01" => -- arithmetic shift right by 1
					t2_s0_mux_out := shift_right(t1(1), 1);
				when "10" => -- arithmetic shift right by 2
					t2_s0_mux_out := shift_right(t1(1), 2);
				when others =>
					t2_s0_mux_out := t1(1); -- default case, no shift
			end case;

			t2(1) := t1(0) - t2_s0_mux_out;

			case i_s1 is
				when '0' => -- mul with k6
					t2_mul1_out := t1(2) * c_k6;
				when '1' => -- mul with k2
					t2_mul1_out := t1(2) * c_k2;
				when others =>
					t2_mul1_out := t1(2); -- default case, no multiplication
			end case;

			case i_s1 is
				when '0' => -- mul with j
					t2(2).re := resize(-t2_mul1_out.im, t2(2).re);
					t2(2).im := t2_mul1_out.re;
				when '1' => -- passthrough
					t2(2) := t2_mul1_out;
				when others =>
					t2(2) := t2_mul1_out; -- default case, passthrough
			end case;

			t2_mul2_out := t1(3) * c_k3;
			t2(3) := t2_mul2_out;

			t2_mul3_out := t1(4) * c_k5;
			t2(4) := t2_mul3_out;

			t2_mul4_out := t1(5) * c_k4;
			t2(5) := t2_mul4_out;

			-- Stage 3
			t3(0) := t2(0);
			t3(5) := t2(1);
			t3(1) := t2(1) + t2(2);
			t3(2) := t2(1) - t2(2);
			t3(3) := t2(3) + t2(5);
			t3(4) := t2(4) + t2(5);

			-- Stage output
			output_X0_wide <= t3(0);

			case i_s0 is
				when "00" =>
					output_X1_wide <= t3(5);
				when "01" =>
					output_X1_wide <= t3(1);
				when others =>
					output_X1_wide <= t3(1)+t3(3);
			end case;

			case i_s1 is
				when '0' =>
					output_X2_wide <= t3(2);
				when '1' =>
					output_X2_wide <= t3(2)+t3(4);
				when others =>
					output_X2_wide <= t3(2); -- default case, passthrough
			end case;

			output_X4_wide <= t3(1) - t3(3);
			output_X3_wide <= t3(2) - t3(4);

		end process PROC_CALC_235;
	end generate GEN_235;

	GEN_23: if G_CAPABILITY = 1 generate
		PROC_CALC_23: process(input_x0_wide, input_x1_wide, input_x2_wide, i_s0)
			variable t1_s0_mux_out, t1_mul1_out : t_cmplx_wide;
			variable t0, t1, t2 : t_cmplx_wide_array;
		begin
			-- Stage 0
			t0(0) := input_x0_wide;
			t0(1) := input_x1_wide + input_x2_wide;
			t0(2) := input_x1_wide - input_x2_wide;

			-- Stage 1
			t1(0) := t0(0) + t0(1);

			case i_s0(0) is
				when '0' => -- no shift
					t1_s0_mux_out := t0(1);
				when '1' => -- arithmetic shift right by 1
					t1_s0_mux_out := shift_right(t0(1), 1);
				when others =>
					t1_s0_mux_out := t1(1); -- default case, no shift
			end case;

			t1(1) := t0(0) - t1_s0_mux_out;

			t1_mul1_out := t0(2) * c_k6;
			-- Apply multiplication by j
			t1(2).re := resize(-t1_mul1_out.im, t1(2).re);
			t1(2).im := t1_mul1_out.re;

			-- Stage 2
			t2(0) := t1(0);
			t2(1) := t1(1) + t1(2);
			t2(2) := t1(1) - t1(2);
			t2(3) := t1(1);

			-- Stage output
			output_X0_wide <= t2(0);
			
			case i_s0(0) is
				when '0' =>
					output_X1_wide <= t2(3);
				when '1' =>
					output_X1_wide <= t2(1);
				when others =>
					output_X1_wide <= t2(3); -- default case
			end case;
			
			output_X2_wide <= t2(2);

		end process PROC_CALC_23;
	end generate GEN_23;

	GEN_2: if G_CAPABILITY = 0 generate
		PROC_CALC_2: process(input_x0_wide, input_x1_wide)
			variable t0, t1 : t_cmplx_wide_array;
		begin
			-- Stage 0
			t0(0) := input_x0_wide;
			t0(1) := input_x1_wide;

			-- Stage 1
			t1(0) := t0(0) + t0(1);
			t1(1) := t0(0) - t0(1);

			-- Stage output
			output_X0_wide <= t1(0);
			output_X1_wide <= t1(1);

		end process PROC_CALC_2;
	end generate GEN_2;

	o_X0 <= resize(output_X0_wide, o_X0);
	o_X1 <= resize(output_X1_wide, o_X1);
	o_X2 <= resize(output_X2_wide, o_X2);
	o_X3 <= resize(output_X3_wide, o_X3);
	o_X4 <= resize(output_X4_wide, o_X4);
	
end architecture rtl;