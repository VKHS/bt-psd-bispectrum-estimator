// Chunk-merge PSD estimator: two memory-based FFT cores plus the segment engine.
//
// Sample stream in, accumulated power spectrum out.  Bit-exact with btpsd.estimator.psd_fixed.
//
//   1. 512 real samples (in_valid / in_data, accepted while in_ready) are loaded into core[cur],
//   2. the core transforms the chunk (scale 1/512),
//   3. from the second chunk on, the segment (older chunk in core[~cur], newer chunk in core[cur]) is merged, windowed,
//      squared and accumulated by psd_engine, then the next chunk is loaded into the core of the older chunk.
//
// After the last chunk has been accepted wait for in_ready (all segments finished), then read the 1024-bin result
// through res_addr / res_data (one cycle latency).  avg_shift = ceil(log2(number of segments)) must be set before the
// last segment finishes.  The accumulators are cleared once after reset (about 1024 cycles).
//
// Reference microarchitecture: loading, transforming and merging are not overlapped; see docs/ARCHITECTURE.md.
`timescale 1ns / 1ps
`include "btpsd_consts.vh"

module bt_psd_top #(
    parameter WL        = 12,
    parameter ROUND     = 1,
    parameter SQ_SERIAL = 1,
    parameter ACCW      = 2 * WL + 8,
    parameter signed [15:0] WM2 = `WIN_HANN_M2,
    parameter signed [15:0] WM1 = `WIN_HANN_M1,
    parameter signed [15:0] W0  = `WIN_HANN_0,
    parameter signed [15:0] WP1 = `WIN_HANN_P1,
    parameter signed [15:0] WP2 = `WIN_HANN_P2
) (
    input                      clk,
    input                      rst_n,
    input      [4:0]           avg_shift,
    input                      in_valid,
    input  signed [WL-1:0]     in_data,
    output                     in_ready,
    output                     busy,
    input      [9:0]           res_addr,
    output     [ACCW-1:0]      res_data
);
    localparam [2:0] T_INIT = 3'd0, T_CLR = 3'd1, T_LOAD = 3'd2, T_FFT = 3'd3, T_SEG = 3'd4;

    reg [2:0] st;
    reg       cur;            // core that receives the chunk being loaded / transformed (the newer chunk of the segment)
    reg       have_prev;
    reg [8:0] ld_cnt;
    reg       c0_start, c1_start, eng_clr, eng_go;

    wire c0_busy, c0_done, c1_busy, c1_done;
    wire [8:0] c0_rd_k, c1_rd_k, eng_a_addr, eng_b_addr;
    wire signed [WL-1:0] c0_re, c0_im, c1_re, c1_im;
    wire eng_busy, eng_done;

    wire ld_now = (st == T_LOAD) && in_valid;
    assign in_ready = (st == T_LOAD);
    assign busy     = (st != T_LOAD);

    fft512_mem #(.WL(WL), .ROUND(ROUND)) u_core0 (
        .clk(clk), .rst_n(rst_n),
        .ld_en(ld_now && !cur), .ld_addr(ld_cnt), .ld_re(in_data), .ld_im({WL{1'b0}}),
        .start(c0_start), .busy(c0_busy), .done(c0_done),
        .rd_k(c0_rd_k), .rd_re(c0_re), .rd_im(c0_im));

    fft512_mem #(.WL(WL), .ROUND(ROUND)) u_core1 (
        .clk(clk), .rst_n(rst_n),
        .ld_en(ld_now && cur), .ld_addr(ld_cnt), .ld_re(in_data), .ld_im({WL{1'b0}}),
        .start(c1_start), .busy(c1_busy), .done(c1_done),
        .rd_k(c1_rd_k), .rd_re(c1_re), .rd_im(c1_im));

    // the engine reads A from the older core and B from the newer core
    assign c0_rd_k = cur ? eng_a_addr : eng_b_addr;
    assign c1_rd_k = cur ? eng_b_addr : eng_a_addr;
    wire signed [WL-1:0] a_re = cur ? c0_re : c1_re;
    wire signed [WL-1:0] a_im = cur ? c0_im : c1_im;
    wire signed [WL-1:0] b_re = cur ? c1_re : c0_re;
    wire signed [WL-1:0] b_im = cur ? c1_im : c0_im;

    psd_engine #(
        .WL(WL), .ROUND(ROUND), .ACCW(ACCW), .SQ_SERIAL(SQ_SERIAL),
        .WM2(WM2), .WM1(WM1), .W0(W0), .WP1(WP1), .WP2(WP2)
    ) u_eng (
        .clk(clk), .rst_n(rst_n), .clr(eng_clr), .go(eng_go), .busy(eng_busy), .done(eng_done),
        .a_addr(eng_a_addr), .a_re(a_re), .a_im(a_im),
        .b_addr(eng_b_addr), .b_re(b_re), .b_im(b_im),
        .avg_shift(avg_shift), .res_addr(res_addr), .res_data(res_data));

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            st <= T_INIT;  cur <= 1'b0;  have_prev <= 1'b0;  ld_cnt <= 9'd0;
            c0_start <= 1'b0;  c1_start <= 1'b0;  eng_clr <= 1'b0;  eng_go <= 1'b0;
        end else begin
            c0_start <= 1'b0;  c1_start <= 1'b0;  eng_clr <= 1'b0;  eng_go <= 1'b0;
            case (st)
                T_INIT: begin
                    eng_clr <= 1'b1;
                    st <= T_CLR;
                end
                T_CLR: begin
                    if (eng_done) st <= T_LOAD;
                end
                T_LOAD: begin
                    if (in_valid) begin
                        if (ld_cnt == 9'd511) begin
                            ld_cnt <= 9'd0;
                            if (cur) c1_start <= 1'b1; else c0_start <= 1'b1;
                            st <= T_FFT;
                        end else begin
                            ld_cnt <= ld_cnt + 9'd1;
                        end
                    end
                end
                T_FFT: begin
                    if (cur ? c1_done : c0_done) begin
                        if (have_prev) begin
                            eng_go <= 1'b1;
                            st <= T_SEG;
                        end else begin
                            have_prev <= 1'b1;
                            cur <= ~cur;
                            st <= T_LOAD;
                        end
                    end
                end
                default: begin   // T_SEG
                    if (eng_done) begin
                        cur <= ~cur;
                        st <= T_LOAD;
                    end
                end
            endcase
        end
    end
endmodule
