// Memory-based 512-point FFT: in-place radix-8 DIF on eight 64-word banks, one butterfly issued per cycle.
//
// Mirrors btpsd.schedule / btpsd.r8fft (bit-exact, scaled by 1/512):
//   address a = 64 a2 + 8 a1 + a0,  bank = (a2 + a1 + a0) mod 8,  row = 8 a2 + a1
//   pass 0 / 1 / 2 combine the eight words that differ in digit a2 / a1 / a0; they always sit in eight different banks.
//   butterfly b = 8 hi + lo of pass p, slot m:  bank = (m + hi + lo) mod 8,
//       row = {m, hi} (pass 0), {hi, m} (pass 1), {hi, lo} (pass 2)
//   twiddle exponent of output k (units of pi/512): pass 0: 2 b k, pass 1: 16 lo k, pass 2: none
//   after pass 2 the word at address (a2, a1, a0) is the spectrum value X[a2 + 8 a1 + 64 a0].
//
// Interface (all synchronous to clk):
//   load    ld_en, ld_addr (natural order), ld_re, ld_im   - only while idle
//   run     start (pulse while idle) -> busy ... done (one-cycle pulse, results are in the banks)
//   unload  rd_k (natural frequency index) -> rd_re / rd_im one cycle later, only while idle
//
// Schedule: LAT = LAT_BFLY + 2 cycles from issue to the cycle in which the written word can be read again.  Pass 1 can
// start LAT - 8 cycles late at most (the producers of its first butterfly are pass-0 butterflies 0, 8, ..., 56), so
// STALL = max(0, LAT - 8) idle cycles are inserted after pass 0; there is no other hazard (btpsd.schedule.schedule_cycles).
`timescale 1ns / 1ps
`include "btpsd_params.vh"

module fft512_mem #(
    parameter WL    = 16,
    parameter ROUND = 1
) (
    input                      clk,
    input                      rst_n,
    // load (natural order)
    input                      ld_en,
    input      [8:0]           ld_addr,
    input  signed [WL-1:0]     ld_re,
    input  signed [WL-1:0]     ld_im,
    // control
    input                      start,
    output reg                 busy,
    output reg                 done,
    // unload (natural frequency order, data one cycle after rd_k)
    input      [8:0]           rd_k,
    output signed [WL-1:0]     rd_re,
    output signed [WL-1:0]     rd_im
);
    localparam LB    = `LAT_BFLY;
    localparam LAT   = LB + 2;
    localparam STALL = (LAT > 8) ? (LAT - 8) : 0;
    localparam DW    = 2 * WL;

    localparam [1:0] ST_IDLE = 2'd0, ST_RUN = 2'd1, ST_DRAIN = 2'd2;

    reg [1:0] st;
    reg       iss_v;                  // a butterfly is issued in this cycle
    reg [1:0] iss_p;
    reg [5:0] iss_b;
    reg [3:0] stall_cnt;
    reg [4:0] drain_cnt;

    wire run_mode = (st != ST_IDLE);

    // ---- issue controller ---------------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            st <= ST_IDLE;  busy <= 1'b0;  done <= 1'b0;
            iss_v <= 1'b0;  iss_p <= 2'd0;  iss_b <= 6'd0;
            stall_cnt <= 4'd0;  drain_cnt <= 5'd0;
        end else begin
            done <= 1'b0;
            case (st)
                ST_IDLE: begin
                    if (start) begin
                        st <= ST_RUN;  busy <= 1'b1;
                        iss_v <= 1'b1;  iss_p <= 2'd0;  iss_b <= 6'd0;
                    end
                end
                ST_RUN: begin
                    if (stall_cnt != 4'd0) begin
                        stall_cnt <= stall_cnt - 4'd1;
                        if (stall_cnt == 4'd1) iss_v <= 1'b1;
                    end else if (iss_b == 6'd63) begin
                        if (iss_p == 2'd2) begin
                            iss_v <= 1'b0;
                            st <= ST_DRAIN;
                            drain_cnt <= LAT - 2;
                        end else begin
                            iss_p <= iss_p + 2'd1;
                            iss_b <= 6'd0;
                            if (iss_p == 2'd0 && STALL != 0) begin
                                iss_v <= 1'b0;
                                stall_cnt <= STALL;
                            end
                        end
                    end else begin
                        iss_b <= iss_b + 6'd1;
                    end
                end
                default: begin   // ST_DRAIN
                    if (drain_cnt == 5'd0) begin
                        st <= ST_IDLE;  busy <= 1'b0;  done <= 1'b1;
                    end else begin
                        drain_cnt <= drain_cnt - 5'd1;
                    end
                end
            endcase
        end
    end

    // ---- address functions -------------------------------------------------------------
    function [5:0] row_fn;
        input [1:0] p;
        input [2:0] m;
        input [2:0] hi;
        input [2:0] lo;
        begin
            case (p)
                2'd0:    row_fn = {m, hi};
                2'd1:    row_fn = {hi, m};
                default: row_fn = {hi, lo};
            endcase
        end
    endfunction

    function [9:0] twq;                      // twiddle exponent, units of pi/512
        input [1:0] p;
        input [5:0] b;
        input [2:0] k;
        reg   [9:0] t;
        begin
            case (p)
                2'd0:    t = ({4'd0, b} * {7'd0, k}) << 1;
                2'd1:    t = ({7'd0, b[2:0]} * {7'd0, k}) << 4;
                default: t = 10'd0;
            endcase
            twq = t;
        end
    endfunction

    // ---- control pipeline ---------------------------------------------------------------
    wire        v1;
    wire [1:0]  p1;
    wire [5:0]  b1;
    delay_line #(.W(9), .N(1)) u_meta1 (.clk(clk), .d({iss_v, iss_p, iss_b}), .q({v1, p1, b1}));

    wire        vw;
    wire [1:0]  pw;
    wire [5:0]  bw;
    delay_line #(.W(9), .N(LB)) u_meta2 (.clk(clk), .d({v1, p1, b1}), .q({vw, pw, bw}));

    // ---- load / unload address decode ------------------------------------------------------
    wire [2:0] ld_bank = ld_addr[8:6] + ld_addr[5:3] + ld_addr[2:0];
    wire [5:0] ld_row  = {ld_addr[8:6], ld_addr[5:3]};

    wire [2:0] ul_bank = rd_k[2:0] + rd_k[5:3] + rd_k[8:6];     // a2 = k[2:0], a1 = k[5:3], a0 = k[8:6]
    wire [5:0] ul_row  = {rd_k[2:0], rd_k[5:3]};
    reg  [2:0] ul_bank_q;
    always @(posedge clk) ul_bank_q <= ul_bank;

    // ---- banks ------------------------------------------------------------------------------
    wire [8*DW-1:0] rq_all;                  // bank j at [j*DW +: DW] = {re, im}
    wire [8*WL-1:0] bf_yr, bf_yi;

    wire [2:0] s_rd = iss_b[5:3] + iss_b[2:0];
    wire [2:0] s_wr = bw[5:3] + bw[2:0];

    genvar j;
    generate
        for (j = 0; j < 8; j = j + 1) begin : g_bank
            localparam [2:0] JB = j;
            wire [2:0]    m_r    = JB - s_rd;
            wire [5:0]    raddr  = run_mode ? row_fn(iss_p, m_r, iss_b[5:3], iss_b[2:0]) : ul_row;
            wire [2:0]    m_w    = JB - s_wr;
            wire [5:0]    wrow   = row_fn(pw, m_w, bw[5:3], bw[2:0]);
            wire [WL-1:0] wd_re  = bf_yr[m_w*WL +: WL];
            wire [WL-1:0] wd_im  = bf_yi[m_w*WL +: WL];
            wire          we     = run_mode ? vw : (ld_en && (ld_bank == JB));
            wire [5:0]    waddr  = run_mode ? wrow : ld_row;
            wire [DW-1:0] wdata  = run_mode ? {wd_re, wd_im} : {ld_re, ld_im};
            ram_1r1w #(.DW(DW), .AW(6)) u_ram (
                .clk(clk), .we(we), .waddr(waddr), .wdata(wdata),
                .raddr(raddr), .rdata(rq_all[j*DW +: DW]));
        end
    endgenerate

    // ---- butterfly input: slot m comes from bank (m + s1) mod 8 -----------------------------------
    wire [2:0] s1 = b1[5:3] + b1[2:0];
    wire [8*WL-1:0] bf_xr, bf_xi;
    wire [69:0]     bf_tq;
    genvar m;
    generate
        for (m = 0; m < 8; m = m + 1) begin : g_slot
            localparam [2:0] MB = m;
            wire [2:0]    bk = MB + s1;
            wire [DW-1:0] w  = rq_all[bk*DW +: DW];
            assign bf_xr[m*WL +: WL] = w[DW-1:WL];
            assign bf_xi[m*WL +: WL] = w[WL-1:0];
            if (m > 0) begin : g_tq
                assign bf_tq[(m-1)*10 +: 10] = twq(p1, b1, MB);
            end
        end
    endgenerate

    r8_butterfly #(.WL(WL), .ROUND(ROUND)) u_bf (
        .clk(clk), .xr(bf_xr), .xi(bf_xi), .tq(bf_tq), .yr(bf_yr), .yi(bf_yi));

    // ---- unload data --------------------------------------------------------------------------------
    wire [DW-1:0] rd_w = rq_all[ul_bank_q*DW +: DW];
    assign rd_re = rd_w[DW-1:WL];
    assign rd_im = rd_w[WL-1:0];
endmodule
