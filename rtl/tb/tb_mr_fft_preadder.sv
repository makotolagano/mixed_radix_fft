module tb_mr_fft_preadder import mr_fft_pkg::*; (
);
    localparam string VCD_FILE = "build/tb_mr_fft_preadder.vcd";

    mr_fft_pkg::cmplx_type x0, x1, x2, x3, x4;
    mr_fft_pkg::cmplx_type X0, X1, X2, X3, X4;

    logic [1:0] s0;
    logic       s1;

    mr_fft_preadder #(
        .cmplx_t(mr_fft_pkg::cmplx_type)
    ) preadder_inst (
        .i_x0(x0),
        .i_x1(x1),
        .i_x2(x2),
        .i_x3(x3),
        .i_x4(x4),
        .i_s0(s0),
        .i_s1(s1),
        .o_X0(X0),
        .o_X1(X1),
        .o_X2(X2),
        .o_X3(X3),
        .o_X4(X4)
    );

    initial begin
        $dumpfile(VCD_FILE);
        $dumpvars(0, tb_mr_fft_preadder);

        s0 = 2'b00;
        s1 = 1'b0;
        x2 = '0;
        x3 = '0;
        x4 = '0;

        // Test case 1
        x0.re = 10; x0.im = 20;
        x1.re = 5;  x1.im = 15;
        #10; // Wait for the sum to be computed
        $display("Test Case 1: x0=(%0d, %0d), x1=(%0d, %0d) => X0=(%0d, %0d), X1=(%0d, %0d)", 
                 x0.re, x0.im, x1.re, x1.im, X0.re, X0.im, X1.re, X1.im);
        
        // Test case 2
        x0.re = -10; x0.im = -20;
        x1.re = -5;  x1.im = -15;
        #10; // Wait for the sum to be computed
        $display("Test Case 2: x0=(%0d, %0d), x1=(%0d, %0d) => X0=(%0d, %0d), X1=(%0d, %0d)", 
                 x0.re, x0.im, x1.re, x1.im, X0.re, X0.im, X1.re, X1.im);
        
        // Test case 3
        x0.re = 32767; x0.im = 32767;
        x1.re = 1;     x1.im = 1;
        #10; // Wait for the sum to be computed
        $display("Test Case 3: x0=(%0d, %0d), x1=(%0d, %0d) => X0=(%0d, %0d), X1=(%0d, %0d)", 
                 x0.re, x0.im, x1.re, x1.im, X0.re, X0.im, X1.re, X1.im);
        
        $finish;
    end

    
endmodule
