// Direct bispectrum estimator: framing with hop 256, one memory-based FFT core and the bispectrum engine.
//
// Sample stream in, accumulated bispectrum B[f1][f2] (f1, f2 < 2**LOGB) out.  Bit-exact with
// btpsd.bispectrum.bispectrum_fixed for 512-sample frames and hop 256 (50 % overlap).
//
//   frame 0   512 samples from the stream into the core; the second half is also kept in a 256-word buffer
//   frame n   the first half is copied from that buffer, 256 new samples complete the frame
//   then      FFT (scale 1/512), bispectrum_engine (window, products, accumulation), next frame
//
// After the last frame has been accepted wait for in_ready, then read the result with res_idx = {f1, f2}
// (one cycle latency).  avg_shift = ceil(log2(number of frames)) must be set before the last frame finishes.
`timescale 1ns / 1ps
`include "btpsd_consts.vh"

module bt_bispectrum_top #(
    parameter WL    = 12,
    parameter ROUND = 1,
    parameter LOGB  = 5,
    parameter BACCW = 2 * WL + 10,
    parameter signed [15:0] WM2 = `WIN_HANN_M2,
    parameter signed [15:0] WM1 = `WIN_HANN_M1,
    parameter signed [15:0] W0  = `WIN_HANN_0,
    parameter signed [15:0] WP1 = `WIN_HANN_P1,
    parameter signed [15:0] WP2 = `WIN_HANN_P2
) (
    input                          clk,
    input                          rst_n,
    input      [4:0]               avg_shift,
    input                          in_valid,
    input  signed [WL-1:0]         in_data,
    output                         in_ready,
    output                         busy,
    input      [2*LOGB-1:0]        res_idx,
    output signed [BACCW-1:0]      res_re,
    output signed [BACCW-1:0]      res_im
);
    localparam [2:0] B_INIT = 3'd0, B_CLR = 3'd1, B_COPY = 3'd2, B_EXT = 3'd3, B_FFT = 3'd4, B_ENG = 3'd5;

    reg [2:0] st;
    reg       first;
    reg [8:0] ld_cnt, cp_cnt;
    reg       cp_v;
    reg [7:0] cp_addr;
    reg       core_start, eng_clr, eng_go;

    wire core_busy, core_done, eng_busy, eng_done;
    wire [8:0] eng_f_addr;
    wire signed [WL-1:0] c_re, c_im;
    wire [WL-1:0] half_q;

    assign in_ready = (st == B_EXT);
    assign busy     = (st != B_EXT);

    wire ext_wr  = (st == B_EXT) && in_valid;
    wire core_ld = cp_v | ext_wr;

    ram_1r1w #(.DW(WL), .AW(8)) u_half (
        .clk(clk), .we(ext_wr && ld_cnt[8]), .waddr(ld_cnt[7:0]), .wdata(in_data),
        .raddr(cp_cnt[7:0]), .rdata(half_q));

    fft512_mem #(.WL(WL), .ROUND(ROUND)) u_core (
        .clk(clk), .rst_n(rst_n),
        .ld_en(core_ld), .ld_addr(cp_v ? {1'b0, cp_addr} : ld_cnt),
        .ld_re(cp_v ? half_q : in_data), .ld_im({WL{1'b0}}),
        .start(core_start), .busy(core_busy), .done(core_done),
        .rd_k(eng_f_addr), .rd_re(c_re), .rd_im(c_im));

    bispectrum_engine #(
        .WL(WL), .ROUND(ROUND), .LOGB(LOGB), .BACCW(BACCW),
        .WM2(WM2), .WM1(WM1), .W0(W0), .WP1(WP1), .WP2(WP2)
    ) u_eng (
        .clk(clk), .rst_n(rst_n), .clr(eng_clr), .go(eng_go), .busy(eng_busy), .done(eng_done),
        .f_addr(eng_f_addr), .f_re(c_re), .f_im(c_im),
        .avg_shift(avg_shift), .res_idx(res_idx), .res_re(res_re), .res_im(res_im));

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            st <= B_INIT;  first <= 1'b1;  ld_cnt <= 9'd0;  cp_cnt <= 9'd0;  cp_v <= 1'b0;  cp_addr <= 8'd0;
            core_start <= 1'b0;  eng_clr <= 1'b0;  eng_go <= 1'b0;
        end else begin
            core_start <= 1'b0;  eng_clr <= 1'b0;  eng_go <= 1'b0;
            cp_v <= 1'b0;
            case (st)
                B_INIT: begin
                    eng_clr <= 1'b1;
                    st <= B_CLR;
                end
                B_CLR: begin
                    if (eng_done) begin
                        ld_cnt <= 9'd0;  cp_cnt <= 9'd0;
                        st <= first ? B_EXT : B_COPY;
                    end
                end
                B_COPY: begin                         // copy the stored second half of the previous frame to addresses 0..255
                    if (cp_cnt < 9'd256) begin
                        cp_v <= 1'b1;  cp_addr <= cp_cnt[7:0];
                        cp_cnt <= cp_cnt + 9'd1;
                    end else begin
                        ld_cnt <= 9'd256;
                        st <= B_EXT;
                    end
                end
                B_EXT: begin
                    if (in_valid) begin
                        if (ld_cnt == 9'd511) begin
                            ld_cnt <= 9'd0;  cp_cnt <= 9'd0;
                            core_start <= 1'b1;
                            st <= B_FFT;
                        end else begin
                            ld_cnt <= ld_cnt + 9'd1;
                        end
                    end
                end
                B_FFT: begin
                    if (core_done) begin
                        eng_go <= 1'b1;
                        st <= B_ENG;
                    end
                end
                default: begin                        // B_ENG
                    if (eng_done) begin
                        first <= 1'b0;
                        st <= B_COPY;
                    end
                end
            endcase
        end
    end
endmodule
