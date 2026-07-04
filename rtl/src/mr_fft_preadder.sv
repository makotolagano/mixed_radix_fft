module mr_fft_preadder import mr_fft_pkg::*; #(
    parameter type cmplx_t = mr_fft_pkg::cmplx_type,
    parameter type cmplx_wide_t = mr_fft_pkg::cmplx_wide_type,
    parameter type cmplx_wide_mult_t = mr_fft_pkg::cmplx_wide_mult_type
) (
    input  cmplx_t i_x0,
    input  cmplx_t i_x1,
    input  cmplx_t i_x2,
    input  cmplx_t i_x3,
    input  cmplx_t i_x4,

    input logic [1:0] i_s0,
    input logic       i_s1,

    output cmplx_t o_X0,
    output cmplx_t o_X1,
    output cmplx_t o_X2,
    output cmplx_t o_X3,
    output cmplx_t o_X4
);

    cmplx_wide_t x0_wide, x1_wide, x2_wide, x3_wide, x4_wide;
    cmplx_wide_t [4:0] t_0;
    cmplx_wide_t [5:0] t_1;
    cmplx_wide_t [5:0] t_2;
    cmplx_wide_t [5:0] t_3;

    cmplx_wide_t t_2_inter, t_out1_r5, t_out1_r3, t_out1_r2, t_out1_r1, t_o_X1, t_o_X2, t_out4, t_out3;
    cmplx_wide_mult_t t_2_mul1, t_2_mul2, t_2_mul3, t_2_mul4;

    // Coefficients
    logic signed [fxp_word_wide_width-1:0] k2, k3, k4, k5, k6;

    assign k2 = mr_fft_pkg::K2;
    assign k3 = mr_fft_pkg::K3;
    assign k4 = mr_fft_pkg::K4;
    assign k5 = mr_fft_pkg::K5;
    assign k6 = mr_fft_pkg::K6;
    

    // Sign-extend and shift left by guard_bits
    always_comb begin : extend_and_shift_inputs
        x0_wide.re = {{guard_bits{i_x0.re[fxp_word_width-1]}}, i_x0.re, {guard_bits{1'b0}}};
        x0_wide.im = {{guard_bits{i_x0.im[fxp_word_width-1]}}, i_x0.im, {guard_bits{1'b0}}};

        x1_wide.re = {{guard_bits{i_x1.re[fxp_word_width-1]}}, i_x1.re, {guard_bits{1'b0}}};
        x1_wide.im = {{guard_bits{i_x1.im[fxp_word_width-1]}}, i_x1.im, {guard_bits{1'b0}}};

        x2_wide.re = {{guard_bits{i_x2.re[fxp_word_width-1]}}, i_x2.re, {guard_bits{1'b0}}};
        x2_wide.im = {{guard_bits{i_x2.im[fxp_word_width-1]}}, i_x2.im, {guard_bits{1'b0}}};

        x3_wide.re = {{guard_bits{i_x3.re[fxp_word_width-1]}}, i_x3.re, {guard_bits{1'b0}}};
        x3_wide.im = {{guard_bits{i_x3.im[fxp_word_width-1]}}, i_x3.im, {guard_bits{1'b0}}};

        x4_wide.re = {{guard_bits{i_x4.re[fxp_word_width-1]}}, i_x4.re, {guard_bits{1'b0}}};
        x4_wide.im = {{guard_bits{i_x4.im[fxp_word_width-1]}}, i_x4.im, {guard_bits{1'b0}}};
    end

    always_comb begin : calc_stages
        // Stage 0
        t_0[0] = x0_wide;

        t_0[1].re = x1_wide.re + x4_wide.re;
        t_0[1].im = x1_wide.im + x4_wide.im;

        t_0[2].re = x2_wide.re + x3_wide.re;
        t_0[2].im = x2_wide.im + x3_wide.im;

        t_0[3].re = x1_wide.re - x4_wide.re;
        t_0[3].im = x1_wide.im - x4_wide.im;

        t_0[4].re = x2_wide.re - x3_wide.re;
        t_0[4].im = x2_wide.im - x3_wide.im;

        // Stage 1
        t_1[0] = t_0[0];

        t_1[1].re = t_0[1].re + t_0[2].re;
        t_1[1].im = t_0[1].im + t_0[2].im;

        t_1[2].re = t_0[1].re - t_0[2].re;
        t_1[2].im = t_0[1].im - t_0[2].im;

        t_1[3] = t_0[3];

        t_1[4] = t_0[4];

        // Stage 2
        t_2[0].re = t_1[0].re + t_1[1].re;
        t_2[0].im = t_1[0].im + t_1[1].im;

        // t_2_inter -> vhdl t2_s0_mux_out
        if (i_s0 == 2'b00) begin // no shift
            t_2_inter = t_1[1];
        end else if (i_s0 == 2'b01) begin // shift > by 1
            t_2_inter.re = t_1[1].re >>> 1; // arithmetic right shift by 1
            t_2_inter.im = t_1[1].im >>> 1;
            // t_2[1].re = {t_2[1].re[fxp_word_wide_width-1], t_2[1].re[fxp_word_wide_width-1:1]};
            // t_2[1].im = {t_2[1].im[fxp_word_wide_width-1], t_2[1].im[fxp_word_wide_width-1:1]};
        end else if (i_s0 == 2'b10) begin // shift > by 2
            t_2_inter.re = t_1[1].re >>> 2; // arithmetic right shift by 2
            t_2_inter.im = t_1[1].im >>> 2;
            // t_2[1].re = {{2{t_2[1].re[fxp_word_wide_width-1]}}, t_2[1].re[fxp_word_wide_width-1:2]};
            // t_2[1].im = {{2{t_2[1].im[fxp_word_wide_width-1]}}, t_2[1].im[fxp_word_wide_width-1:2]};
        end else begin
            t_2_inter = t_1[1];
        end
        t_2[1].re = t_1[0].re - t_2_inter.re;
        t_2[1].im = t_1[0].im - t_2_inter.im;

        if (i_s1 == 1'b1) begin
            t_2_mul1.re = t_1[2].re * k6;
            t_2_mul1.im = t_1[2].im * k6;
        end else begin
            t_2_mul1.re = t_1[2].re * k2;
            t_2_mul1.im = t_1[2].im * k2;
        end
        t_2[2].re = t_2_mul1.re[fxp_word_wide_width*2-fxp_int_wide_width-1 -: fxp_word_wide_width];
        t_2[2].im = t_2_mul1.im[fxp_word_wide_width*2-fxp_int_wide_width-1 -: fxp_word_wide_width];

        if (i_s1 == 1'b0) begin
            t_2[2].re = t_2[2].im;
            t_2[2].im = -t_2[2].re;
        end

        t_2_mul2.re = t_1[3].re * k3;
        t_2_mul2.im = t_1[3].im * k3;
        t_2[3].re = -t_2_mul2.im[fxp_word_wide_width*2-fxp_int_wide_width-1 -: fxp_word_wide_width];
        t_2[3].im = t_2_mul2.re[fxp_word_wide_width*2-fxp_int_wide_width-1 -: fxp_word_wide_width];

        t_2_mul3.re = t_1[4].re * k5;
        t_2_mul3.im = t_1[4].im * k5;
        t_2[4].re = -t_2_mul3.im[fxp_word_wide_width*2-fxp_int_wide_width-1 -: fxp_word_wide_width];
        t_2[4].im = t_2_mul3.re[fxp_word_wide_width*2-fxp_int_wide_width-1 -: fxp_word_wide_width];

        t_2_mul4.re = t_1[4].re * k4;
        t_2_mul4.im = t_1[4].im * k4;
        t_2[5].re = -t_2_mul4.im[fxp_word_wide_width*2-fxp_int_wide_width-1 -: fxp_word_wide_width];
        t_2[5].im = t_2_mul4.re[fxp_word_wide_width*2-fxp_int_wide_width-1 -: fxp_word_wide_width];
    
        // Stage 3
        t_3[0] = t_2[0];
        t_3[5] = t_2[1];
        
        t_3[1].re = t_2[1].re + t_2[2].re;
        t_3[1].im = t_2[1].im + t_2[2].im;

        t_3[2].re = t_2[1].re - t_2[2].re;
        t_3[2].im = t_2[1].im - t_2[2].im;

        t_3[3].re = t_2[3].re + t_2[5].re;
        t_3[3].im = t_2[3].im + t_2[5].im;

        t_3[4].re = t_2[4].re + t_2[5].re;
        t_3[4].im = t_2[4].im + t_2[5].im;

        // Output assignment with rounding and saturation
        o_X0.re = t_3[0].re[fxp_word_wide_width-1 -: fxp_word_width];
        o_X0.im = t_3[0].im[fxp_word_wide_width-1 -: fxp_word_width];

        t_out1_r5.re = t_3[1].re + t_3[3].re;
        t_out1_r5.im = t_3[1].im + t_3[3].im;

        t_out1_r3 = t_3[1];
        t_out1_r2 = t_3[5];

        case (i_s0)
            2'b00: t_o_X1 = t_out1_r2;
            2'b01: t_o_X1 = t_out1_r3;
            default: t_o_X1 = t_out1_r5;
        endcase

        o_X1.re = t_o_X1.re[fxp_word_wide_width-1 -: fxp_word_width];
        o_X1.im = t_o_X1.im[fxp_word_wide_width-1 -: fxp_word_width];
        
        t_out1_r1.re = t_3[2].re + t_3[4].re;
        t_out1_r1.im = t_3[2].im + t_3[4].im;

        if (i_s1 == 1'b0) begin
            t_o_X2.re = t_3[2].im;
            t_o_X2.im = t_3[2].re;
        end else begin
            t_o_X2.re = t_out1_r1.re;
            t_o_X2.im = t_out1_r1.im;
        end
        
        o_X2.re = t_o_X2.re[fxp_word_wide_width-1 -: fxp_word_width];
        o_X2.im = t_o_X2.im[fxp_word_wide_width-1 -: fxp_word_width];

        t_out4.re = t_3[1].re - t_3[3].re;
        t_out4.im = t_3[1].im - t_3[3].im;

        o_X4.re = t_out4.re[fxp_word_wide_width-1 -: fxp_word_width];
        o_X4.im = t_out4.im[fxp_word_wide_width-1 -: fxp_word_width];

        t_out3.re = t_3[2].re - t_3[4].re;
        t_out3.im = t_3[2].im - t_3[4].im;

        o_X3.re = t_out3.re[fxp_word_wide_width-1 -: fxp_word_width];
        o_X3.im = t_out3.im[fxp_word_wide_width-1 -: fxp_word_width];

    end

endmodule
