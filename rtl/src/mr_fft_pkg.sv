package mr_fft_pkg;

    localparam DATA_WIDTH = 16;

    localparam int fxp_word_width = 8;
    localparam int fxp_frac_width = 6;
    localparam int guard_bits = 3;
    localparam int fxp_word_wide_width = fxp_word_width + 2*guard_bits;
    localparam int fxp_frac_wide_width = fxp_frac_width + guard_bits;
    localparam int fxp_int_wide_width  = fxp_word_wide_width - fxp_frac_wide_width;

    typedef struct packed {
        logic signed [fxp_word_width-1:0] re; 
        logic signed [fxp_word_width-1:0] im; 
    } cmplx_type;

    typedef struct packed {
        logic signed [fxp_word_wide_width-1:0] re; 
        logic signed [fxp_word_wide_width-1:0] im; 
    } cmplx_wide_type;

    typedef struct packed {
        logic signed [fxp_word_wide_width*2-1:0] re; 
        logic signed [fxp_word_wide_width*2-1:0] im; 
    } cmplx_wide_mult_type;

    localparam real PI = 3.14159265358979323846;

    function automatic logic signed [fxp_word_wide_width-1:0] fxp_from_real(real value);
        real scaled;
        begin
            scaled = value * (2.0 ** fxp_frac_wide_width);
            fxp_from_real = $rtoi((scaled >= 0.0) ? (scaled + 0.5) : (scaled - 0.5));
        end
    endfunction

    localparam logic signed [fxp_word_wide_width-1:0] K2_RE =
        fxp_from_real(0.5 * ($cos(2.0 * PI / 5.0) - $cos(4.0 * PI / 5.0)));

    localparam logic signed [fxp_word_wide_width-1:0] K3_IM =
        fxp_from_real($sin(4.0 * PI / 5.0) - $sin(2.0 * PI / 5.0));

    localparam logic signed [fxp_word_wide_width-1:0] K4_IM =
        fxp_from_real(-$sin(4.0 * PI / 5.0));

    localparam logic signed [fxp_word_wide_width-1:0] K5_IM =
        fxp_from_real($sin(4.0 * PI / 5.0) + $sin(2.0 * PI / 5.0));

    localparam logic signed [fxp_word_wide_width-1:0] K6_RE =
        fxp_from_real(-$sqrt(3.0) / 2.0);

    localparam signed [fxp_word_wide_width-1:0] K1 = '0;
    localparam signed [fxp_word_wide_width-1:0] K2 = K2_RE;
    localparam signed [fxp_word_wide_width-1:0] K3 = K3_IM;
    localparam signed [fxp_word_wide_width-1:0] K4 = K4_IM;
    localparam signed [fxp_word_wide_width-1:0] K5 = K5_IM;
    localparam signed [fxp_word_wide_width-1:0] K6 = K6_RE;



endpackage
