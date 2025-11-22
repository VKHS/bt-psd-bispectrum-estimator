// bispectrum_engine.v
// Core triple-product for bispectrum:
//   B = F1 * F2 * conj(F3)  (F3 = F(f1+f2)) :contentReference[oaicite:14]{index=14}
//
// Fixed-point complex multipliers from complex_fixed.v.
// This module just does one triple product at a time; you wrap it
// with an accumulator over segments if you want full BT averaging.

`ifndef BISPECTRUM_ENGINE_V
`define BISPECTRUM_ENGINE_V

`include "complex_fixed.v"

module bispectrum_engine #(
    parameter DATA_W = 16,
    parameter FRAC_W = 14
)(
    input  wire                     clk,
    input  wire                     rst_n,
    input  wire                     in_valid,
    input  wire signed [DATA_W-1:0] f1_r,
    input  wire signed [DATA_W-1:0] f1_i,
    input  wire signed [DATA_W-1:0] f2_r,
    input  wire signed [DATA_W-1:0] f2_i,
    input  wire signed [DATA_W-1:0] f3_r,
    input  wire signed [DATA_W-1:0] f3_i,
    output reg                      out_valid,
    output reg  signed [DATA_W-1:0] b_r,
    output reg  signed [DATA_W-1:0] b_i
);
    // Stage 1: t = F1 * F2
    wire signed [DATA_W-1:0] t_r, t_i;

    complex_mul #(
        .DATA_W (DATA_W),
        .FRAC_W (FRAC_W)
    ) mul12 (
        .ar (f1_r), .ai (f1_i),
        .br (f2_r), .bi (f2_i),
        .yr (t_r),  .yi (t_i)
    );

    // Stage 2: B = t * conj(F3)
    wire signed [DATA_W-1:0] b_r_w, b_i_w;

    complex_mul #(
        .DATA_W (DATA_W),
        .FRAC_W (FRAC_W)
    ) mul3 (
        .ar (t_r),     .ai (t_i),
        .br (f3_r),    .bi (-f3_i), // conjugate
        .yr (b_r_w),   .yi (b_i_w)
    );

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            out_valid <= 1'b0;
            b_r <= '0;
            b_i <= '0;
        end else begin
            if (in_valid) begin
                b_r      <= b_r_w;
                b_i      <= b_i_w;
                out_valid <= 1'b1;
            end else begin
                out_valid <= 1'b0;
            end
        end
    end
endmodule

`endif
