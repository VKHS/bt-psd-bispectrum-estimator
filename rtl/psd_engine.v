// PSD segment engine: chunk merge, frequency-domain window, magnitude squared and segment accumulation.
//
// Reference (sequential) microarchitecture of the datapath of paper Figs. 3-5, bit-exact with
// btpsd.estimator.psd_fixed.  One segment is two consecutive 512-sample chunk spectra A (older) and B (newer),
// both held in fft512_mem cores and read through their unload ports (scale 1/512).
//
//   merge   X(2m)   = sat((A(m) + B(m) + 1) >> 1)
//           X(2m+1) = sat(round(-j (-1)^m sum_i h_i (-1)^(m+d_i) D(m+d_i) / 2^(FD_FRAC+1))),   D = A - B,  d_i = -2 .. 3
//           indices wrap modulo 512; the six taps h_i are the integer constants FD_H0..FD_H5 (paper Eq. 6)
//   window  Xw(k) = sat(round(sum_d W(d) X(k-d) / 2^WIN_FRAC)),  d = -2 .. 2, indices wrap modulo 1024
//   square  P(k) = Xw_re^2 + Xw_im^2   (sum-of-odds squarer when SQ_SERIAL = 1)
//   accumulate  ACC(k) += P(k) over segments;  result(k) = round(ACC(k) / 2^avg_shift)
//
// This engine processes one bin at a time to keep the control simple (roughly 3 x 6 cycles per merge step and
// 3 x 5 cycles plus the squarer time per window bin).  It is the functional reference for a streaming version;
// the throughput of the paper (Table 9) is not claimed for it.  The FFT core is the high-throughput part.
//
//   clr    pulse while idle: zero the accumulators (1024 cycles)       go    pulse while idle: process one segment
//   busy   high while clearing or processing                           done  one-cycle pulse when finished
//   res_addr -> res_data: shifted accumulator, valid one cycle after res_addr, only while idle
`timescale 1ns / 1ps
`include "btpsd_consts.vh"

module psd_engine #(
    parameter WL        = 12,
    parameter ROUND     = 1,
    parameter ACCW      = 2 * WL + 8,
    parameter SQ_SERIAL = 1,
    parameter signed [15:0] H0  = `FD_H0,
    parameter signed [15:0] H1  = `FD_H1,
    parameter signed [15:0] H2  = `FD_H2,
    parameter signed [15:0] H3  = `FD_H3,
    parameter signed [15:0] H4  = `FD_H4,
    parameter signed [15:0] H5  = `FD_H5,
    parameter signed [15:0] WM2 = `WIN_HANN_M2,
    parameter signed [15:0] WM1 = `WIN_HANN_M1,
    parameter signed [15:0] W0  = `WIN_HANN_0,
    parameter signed [15:0] WP1 = `WIN_HANN_P1,
    parameter signed [15:0] WP2 = `WIN_HANN_P2
) (
    input                      clk,
    input                      rst_n,
    input                      clr,
    input                      go,
    output reg                 busy,
    output reg                 done,
    // spectrum of the older chunk (A) and of the newer chunk (B), natural frequency order, data one cycle after the address
    output reg [8:0]           a_addr,
    input  signed [WL-1:0]     a_re,
    input  signed [WL-1:0]     a_im,
    output reg [8:0]           b_addr,
    input  signed [WL-1:0]     b_re,
    input  signed [WL-1:0]     b_im,
    // result read-out
    input      [4:0]           avg_shift,
    input      [9:0]           res_addr,
    output     [ACCW-1:0]      res_data
);
    localparam FD_SH  = `FD_FRAC + 1;
    localparam WIN_SH = `WIN_FRAC;
    localparam signed [WL+19:0] MAXV = (1 <<< (WL - 1)) - 1;
    localparam signed [WL+19:0] MINV = -(1 <<< (WL - 1));

    localparam [3:0] S_IDLE = 4'd0,  S_CLR = 4'd1,  M_A = 4'd2,  M_B = 4'd3,  M_C = 4'd4,  M_F = 4'd5,
                     W_A = 4'd6,  W_B = 4'd7,  W_C = 4'd8,  W_F = 4'd9,  W_G = 4'd10, W_S = 4'd11,
                     W_ACC = 4'd12, S_END = 4'd13;

    reg [3:0] st;
    reg [8:0] m;
    reg [2:0] t;
    reg [9:0] k;
    reg [9:0] cnt;
    reg [8:0] idx_q;

    // ---- memories -----------------------------------------------------------------------------
    reg                 x_we;
    reg  [9:0]          x_waddr;
    reg  [2*WL-1:0]     x_wdata;
    reg  [9:0]          x_raddr;
    wire [2*WL-1:0]     x_rdata;
    ram_1r1w #(.DW(2*WL), .AW(10)) u_xseg (
        .clk(clk), .we(x_we), .waddr(x_waddr), .wdata(x_wdata), .raddr(x_raddr), .rdata(x_rdata));

    reg                 acc_we;
    reg  [9:0]          acc_waddr;
    reg  [ACCW-1:0]     acc_wdata;
    wire [9:0]          acc_raddr = (st == S_IDLE) ? res_addr : k;
    wire [ACCW-1:0]     acc_rdata;
    ram_1r1w #(.DW(ACCW), .AW(10)) u_acc (
        .clk(clk), .we(acc_we), .waddr(acc_waddr), .wdata(acc_wdata), .raddr(acc_raddr), .rdata(acc_rdata));

    // ---- result scaling ---------------------------------------------------------------------------
    wire [ACCW-1:0] rnd_c = (ROUND && avg_shift != 5'd0) ? ({{(ACCW-1){1'b0}}, 1'b1} << (avg_shift - 5'd1)) : {ACCW{1'b0}};
    assign res_data = (acc_rdata + rnd_c) >> avg_shift;

    // ---- squarer ---------------------------------------------------------------------------------------
    reg                 sq_start;
    reg  signed [WL-1:0] sq_re, sq_im;
    wire                sq_busy, sq_done;
    wire [2*WL-1:0]     sq_p;
    mag_sq_odd #(.WL(WL), .SERIAL(SQ_SERIAL)) u_sq (
        .clk(clk), .rst_n(rst_n), .start(sq_start), .re(sq_re), .im(sq_im),
        .busy(sq_busy), .done(sq_done), .p(sq_p));

    // ---- merge arithmetic -------------------------------------------------------------------------------
    wire [8:0] idx = m + {6'd0, t} - 9'd2;                      // (m + d) mod 512, d = t - 2

    wire signed [WL:0] dre = a_re - b_re;
    wire signed [WL:0] dim = a_im - b_im;
    wire signed [WL:0] dpr = idx_q[0] ? -dre : dre;              // (-1)^idx * D
    wire signed [WL:0] dpi = idx_q[0] ? -dim : dim;

    function signed [15:0] hsel;
        input [2:0] tt;
        case (tt)
            3'd0: hsel = H0;
            3'd1: hsel = H1;
            3'd2: hsel = H2;
            3'd3: hsel = H3;
            3'd4: hsel = H4;
            default: hsel = H5;
        endcase
    endfunction

    function signed [15:0] wsel;
        input [2:0] tt;
        case (tt)
            3'd0: wsel = WM2;
            3'd1: wsel = WM1;
            3'd2: wsel = W0;
            3'd3: wsel = WP1;
            default: wsel = WP2;
        endcase
    endfunction

    wire signed [15:0]    hc  = hsel(t);
    wire signed [WL+16:0] prr = dpr * hc;
    wire signed [WL+16:0] pri = dpi * hc;

    wire signed [WL+1:0]  esr = a_re + b_re;
    wire signed [WL+1:0]  esi = a_im + b_im;
    wire signed [WL+1:0]  ehr = ROUND ? (esr + 1) : esr;
    wire signed [WL+1:0]  ehi = ROUND ? (esi + 1) : esi;
    wire signed [WL+1:0]  er2 = ehr >>> 1;
    wire signed [WL+1:0]  ei2 = ehi >>> 1;

    reg signed [WL+19:0]  facc_r, facc_i;
    wire signed [WL+19:0] ow_r = m[0] ? -facc_i : facc_i;       // -j (-1)^m acc = ((-1)^m acc_i, -(-1)^m acc_r)
    wire signed [WL+19:0] ow_i = m[0] ? facc_r : -facc_r;

    // ---- window arithmetic --------------------------------------------------------------------------------
    wire [9:0] idxw = k + 10'd2 - {7'd0, t};                    // (k - d) mod 1024, d = t - 2
    wire signed [WL-1:0]  xr_w = x_rdata[2*WL-1:WL];
    wire signed [WL-1:0]  xi_w = x_rdata[WL-1:0];
    wire signed [15:0]    wc   = wsel(t);
    wire signed [WL+15:0] pwr  = xr_w * wc;
    wire signed [WL+15:0] pwi  = xi_w * wc;
    reg signed [WL+19:0]  wacc_r, wacc_i;

    // round (half up), shift right by sh, saturate to WL bits
    function signed [WL-1:0] rsat;
        input signed [WL+19:0] v;
        input integer          sh;
        reg signed [WL+19:0]   tv;
        begin
            tv = ROUND ? (v + (1 <<< (sh - 1))) : v;
            tv = tv >>> sh;
            if (tv > MAXV)      rsat = MAXV[WL-1:0];
            else if (tv < MINV) rsat = MINV[WL-1:0];
            else                rsat = tv[WL-1:0];
        end
    endfunction

    // ---- control ---------------------------------------------------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            st <= S_IDLE;  busy <= 1'b0;  done <= 1'b0;
            m <= 0;  t <= 0;  k <= 0;  cnt <= 0;  idx_q <= 0;
            a_addr <= 0;  b_addr <= 0;
            x_we <= 0;  x_waddr <= 0;  x_wdata <= 0;  x_raddr <= 0;
            acc_we <= 0;  acc_waddr <= 0;  acc_wdata <= 0;
            sq_start <= 0;  sq_re <= 0;  sq_im <= 0;
            facc_r <= 0;  facc_i <= 0;  wacc_r <= 0;  wacc_i <= 0;
        end else begin
            done     <= 1'b0;
            x_we     <= 1'b0;
            acc_we   <= 1'b0;
            sq_start <= 1'b0;
            case (st)
                S_IDLE: begin
                    if (clr) begin
                        st <= S_CLR;  busy <= 1'b1;  cnt <= 10'd0;
                    end else if (go) begin
                        st <= M_A;  busy <= 1'b1;
                        m <= 9'd0;  t <= 3'd0;
                        facc_r <= 0;  facc_i <= 0;
                    end
                end
                S_CLR: begin
                    acc_we <= 1'b1;  acc_waddr <= cnt;  acc_wdata <= {ACCW{1'b0}};
                    cnt <= cnt + 10'd1;
                    if (cnt == 10'd1023) st <= S_END;
                end
                // ---- merge: six taps per m -------------------------------------------------
                M_A: begin
                    a_addr <= idx;  b_addr <= idx;  idx_q <= idx;
                    st <= M_B;
                end
                M_B: st <= M_C;
                M_C: begin
                    facc_r <= facc_r + prr;
                    facc_i <= facc_i + pri;
                    if (t == 3'd2) begin                      // d = 0: A(m), B(m) are the even-bin operands
                        x_we <= 1'b1;  x_waddr <= {m, 1'b0};
                        x_wdata <= {er2[WL-1:0], ei2[WL-1:0]};
                    end
                    if (t == 3'd5) st <= M_F;
                    else begin t <= t + 3'd1;  st <= M_A; end
                end
                M_F: begin
                    x_we <= 1'b1;  x_waddr <= {m, 1'b1};
                    x_wdata <= {rsat(ow_r, FD_SH), rsat(ow_i, FD_SH)};
                    facc_r <= 0;  facc_i <= 0;  t <= 3'd0;
                    if (m == 9'd511) begin
                        k <= 10'd0;  wacc_r <= 0;  wacc_i <= 0;
                        st <= W_A;
                    end else begin
                        m <= m + 9'd1;
                        st <= M_A;
                    end
                end
                // ---- window, square, accumulate: five taps per bin -------------------------------------------
                W_A: begin
                    x_raddr <= idxw;
                    st <= W_B;
                end
                W_B: st <= W_C;
                W_C: begin
                    wacc_r <= wacc_r + pwr;
                    wacc_i <= wacc_i + pwi;
                    if (t == 3'd4) st <= W_F;
                    else begin t <= t + 3'd1;  st <= W_A; end
                end
                W_F: begin
                    sq_re <= rsat(wacc_r, WIN_SH);
                    sq_im <= rsat(wacc_i, WIN_SH);
                    wacc_r <= 0;  wacc_i <= 0;  t <= 3'd0;
                    st <= W_G;
                end
                W_G: begin
                    sq_start <= 1'b1;
                    st <= W_S;
                end
                W_S: begin
                    if (sq_done) st <= W_ACC;
                end
                W_ACC: begin
                    acc_we <= 1'b1;  acc_waddr <= k;  acc_wdata <= acc_rdata + {{(ACCW-2*WL){1'b0}}, sq_p};
                    if (k == 10'd1023) st <= S_END;
                    else begin k <= k + 10'd1;  st <= W_A; end
                end
                default: begin                                  // S_END
                    busy <= 1'b0;  done <= 1'b1;  st <= S_IDLE;
                end
            endcase
        end
    end
endmodule
