// Magnitude squared  p = re^2 + im^2  by the sum-of-odd-numbers squarer of the paper (SERIAL = 1) or by multipliers.
//
// n^2 = 1 + 3 + 5 + ... + (2|n| - 1): one adder and a down counter per operand, no multiplier.  The real and the
// imaginary squarer run in parallel, so a serial evaluation takes max(|re|, |im|) cycles (btpsd.squarer.mag_sq, units = 2)
// plus two control cycles.  With SERIAL = 0 the result is registered after one cycle.
//
//   start  pulse while !busy, with re / im stable in that cycle
//   done   one-cycle pulse, p valid from the cycle in which done is high until the next start
`timescale 1ns / 1ps

module mag_sq_odd #(
    parameter WL     = 16,
    parameter SERIAL = 1
) (
    input                    clk,
    input                    rst_n,
    input                    start,
    input  signed [WL-1:0]   re,
    input  signed [WL-1:0]   im,
    output reg               busy,
    output reg               done,
    output reg [2*WL-1:0]    p
);
    generate
        if (SERIAL) begin : g_serial
            reg [WL-1:0]   cnt_r, cnt_i;          // |re|, |im| remaining
            reg [2*WL-1:0] odd_r, odd_i;          // next odd number
            reg [2*WL-1:0] acc_r, acc_i;

            wire [WL-1:0] abs_r = re[WL-1] ? (~re + 1'b1) : re;    // |-2^(WL-1)| = 2^(WL-1) fits as unsigned
            wire [WL-1:0] abs_i = im[WL-1] ? (~im + 1'b1) : im;

            always @(posedge clk or negedge rst_n) begin
                if (!rst_n) begin
                    busy <= 1'b0;  done <= 1'b0;  p <= {2*WL{1'b0}};
                    cnt_r <= 0;  cnt_i <= 0;  odd_r <= 0;  odd_i <= 0;  acc_r <= 0;  acc_i <= 0;
                end else begin
                    done <= 1'b0;
                    if (!busy) begin
                        if (start) begin
                            busy  <= 1'b1;
                            cnt_r <= abs_r;  cnt_i <= abs_i;
                            odd_r <= 1;      odd_i <= 1;
                            acc_r <= 0;      acc_i <= 0;
                        end
                    end else begin
                        if (cnt_r == 0 && cnt_i == 0) begin
                            p <= acc_r + acc_i;
                            done <= 1'b1;
                            busy <= 1'b0;
                        end else begin
                            if (cnt_r != 0) begin
                                acc_r <= acc_r + odd_r;  odd_r <= odd_r + 2;  cnt_r <= cnt_r - 1'b1;
                            end
                            if (cnt_i != 0) begin
                                acc_i <= acc_i + odd_i;  odd_i <= odd_i + 2;  cnt_i <= cnt_i - 1'b1;
                            end
                        end
                    end
                end
            end
        end else begin : g_mult
            always @(posedge clk or negedge rst_n) begin
                if (!rst_n) begin
                    busy <= 1'b0;  done <= 1'b0;  p <= {2*WL{1'b0}};
                end else begin
                    done <= start;
                    busy <= 1'b0;
                    if (start) p <= re * re + im * im;
                end
            end
        end
    endgenerate
endmodule
