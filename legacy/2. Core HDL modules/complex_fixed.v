// complex_fixed.v
// Basic fixed-point complex arithmetic blocks.
// Parameterized signed Q(FRAC_W) format.
//
// Re/Im are signed DATA_W-wide values, interpreted as fixed-point
// with FRAC_W fractional bits. No explicit saturation here – for
// real hardware you may add saturation logic or rely on safe scaling.

`ifndef COMPLEX_FIXED_V
`define COMPLEX_FIXED_V

module complex_add #(
    parameter DATA_W = 16
)(
    input  wire signed [DATA_W-1:0] ar,
    input  wire signed [DATA_W-1:0] ai,
    input  wire signed [DATA_W-1:0] br,
    input  wire signed [DATA_W-1:0] bi,
    output wire signed [DATA_W-1:0] yr,
    output wire signed [DATA_W-1:0] yi
);
    assign yr = ar + br;
    assign yi = ai + bi;
endmodule


module complex_sub #(
    parameter DATA_W = 16
)(
    input  wire signed [DATA_W-1:0] ar,
    input  wire signed [DATA_W-1:0] ai,
    input  wire signed [DATA_W-1:0] br,
    input  wire signed [DATA_W-1:0] bi,
    output wire signed [DATA_W-1:0] yr,
    output wire signed [DATA_W-1:0] yi
);
    assign yr = ar - br;
    assign yi = ai - bi;
endmodule


module complex_mul #(
    parameter DATA_W = 16,
    parameter FRAC_W = 14   // number of fractional bits
)(
    input  wire signed [DATA_W-1:0] ar,
    input  wire signed [DATA_W-1:0] ai,
    input  wire signed [DATA_W-1:0] br,
    input  wire signed [DATA_W-1:0] bi,
    output wire signed [DATA_W-1:0] yr,
    output wire signed [DATA_W-1:0] yi
);
    // Intermediate full-precision products (2*DATA_W bits)
    wire signed [2*DATA_W-1:0] p1 = ar * br;
    wire signed [2*DATA_W-1:0] p2 = ai * bi;
    wire signed [2*DATA_W-1:0] p3 = ar * bi;
    wire signed [2*DATA_W-1:0] p4 = ai * br;

    // Real = (ar*br - ai*bi) >> FRAC_W
    // Imag = (ar*bi + ai*br) >> FRAC_W
    wire signed [2*DATA_W:0] real_full = p1 - p2;
    wire signed [2*DATA_W:0] imag_full = p3 + p4;

    assign yr = real_full >>> FRAC_W;
    assign yi = imag_full >>> FRAC_W;
endmodule

`endif
