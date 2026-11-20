// Direct bispectrum engine: frequency-domain window, then B(f1, f2) += P12 * conj(F(f1 + f2)) over frames.
//
// Reference (sequential) microarchitecture, bit-exact with btpsd.bispectrum.bispectrum_fixed.  One 512-sample frame
// spectrum F (scale 1/512) is read from an fft512_mem unload port.
//
//   window    Fw(k) = sat(round(sum_d W(d) F(k-d) / 2^WIN_FRAC)),  d = -2..2, indices modulo 512, for k = 0 .. 2*MAXB-1
//   product   P12 = round(Fw(f1) * Fw(f2) / 2^(WL-1))   (complex, real and imaginary part rounded separately)
//   triple    T = P12 * conj(Fw((f1 + f2) mod 512))      (exact)
//   accumulate  ACC[f1][f2] += T   over frames, f1, f2 = 0 .. MAXB-1;     result = round(ACC / 2^avg_shift)
//
// The windowed spectrum is kept in two identical RAMs so that Fw(f2) and Fw(f1 + f2) are read in the same cycle.
// One (f1, f2) pair takes four cycles; MAXB = 2**LOGB.  The long multiply-accumulate chains are not pipelined
// (functional reference, see docs/ARCHITECTURE.md).
//
//   clr    pulse while idle: zero the accumulators      go    pulse while idle: process one frame
//   busy / done as in psd_engine                        res_idx -> res_re / res_im one cycle later, only while idle
`timescale 1ns / 1ps
`include "btpsd_consts.vh"

module bispectrum_engine #(
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
    input                          clr,
    input                          go,
    output reg                     busy,
    output reg                     done,
    output reg [8:0]               f_addr,
    input  signed [WL-1:0]         f_re,
    input  signed [WL-1:0]         f_im,
    input      [4:0]               avg_shift,
    input      [2*LOGB-1:0]        res_idx,
    output signed [BACCW-1:0]      res_re,
    output signed [BACCW-1:0]      res_im
);
    localparam MAXB  = 1 << LOGB;
    localparam NBINS = 2 * MAXB;
    localparam SHP   = WL - 1;
    localparam signed [WL+19:0] MAXV = (1 <<< (WL - 1)) - 1;
    localparam signed [WL+19:0] MINV = -(1 <<< (WL - 1));

    localparam [3:0] S_IDLE = 4'd0, S_CLR = 4'd1, W_A = 4'd2, W_B = 4'd3, W_C = 4'd4, W_F = 4'd5,
                     R_A = 4'd6, R_B = 4'd7, R_C = 4'd8, T_A = 4'd9, T_B = 4'd10, T_C = 4'd11, T_D = 4'd12,
                     S_END = 4'd13;

    reg [3:0]      st;
    reg [2:0]      t;
    reg [8:0]      k;
    reg [LOGB-1:0] f1, f2;
    reg [2*LOGB-1:0] cnt;

    // ---- windowed spectrum: two copies --------------------------------------------------------
    reg            s_we;
    reg  [8:0]     s_waddr;
    reg  [2*WL-1:0] s_wdata;
    reg  [8:0]     ra_addr, rb_addr;
    wire [2*WL-1:0] sa_q, sb_q;
    ram_1r1w #(.DW(2*WL), .AW(9)) u_sa (.clk(clk), .we(s_we), .waddr(s_waddr), .wdata(s_wdata), .raddr(ra_addr), .rdata(sa_q));
    ram_1r1w #(.DW(2*WL), .AW(9)) u_sb (.clk(clk), .we(s_we), .waddr(s_waddr), .wdata(s_wdata), .raddr(rb_addr), .rdata(sb_q));

    // ---- accumulators (real and imaginary) ------------------------------------------------------
    reg                    acc_we;
    reg  [2*LOGB-1:0]      acc_waddr, acc_ra;
    reg  [BACCW-1:0]       accr_wdata, acci_wdata;
    wire [2*LOGB-1:0]      acc_raddr = (st == S_IDLE) ? res_idx : acc_ra;
    wire signed [BACCW-1:0] accr_q, acci_q;
    ram_1r1w #(.DW(BACCW), .AW(2*LOGB)) u_accr (.clk(clk), .we(acc_we), .waddr(acc_waddr), .wdata(accr_wdata), .raddr(acc_raddr), .rdata(accr_q));
    ram_1r1w #(.DW(BACCW), .AW(2*LOGB)) u_acci (.clk(clk), .we(acc_we), .waddr(acc_waddr), .wdata(acci_wdata), .raddr(acc_raddr), .rdata(acci_q));

    wire signed [BACCW-1:0] rnd_c = (ROUND && avg_shift != 5'd0) ? (1 <<< (avg_shift - 5'd1)) : 0;
    assign res_re = (accr_q + rnd_c) >>> avg_shift;
    assign res_im = (acci_q + rnd_c) >>> avg_shift;

    // ---- window arithmetic -----------------------------------------------------------------------
    wire [8:0] idxw = k + 9'd2 - {6'd0, t};                     // (k - d) mod 512, d = t - 2
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
    wire signed [15:0]    wc  = wsel(t);
    wire signed [WL+15:0] pwr = f_re * wc;
    wire signed [WL+15:0] pwi = f_im * wc;
    reg signed [WL+19:0]  wacc_r, wacc_i;

    function signed [WL-1:0] rsat_win;
        input signed [WL+19:0] v;
        reg signed [WL+19:0]   tv;
        begin
            tv = ROUND ? (v + (1 <<< (`WIN_FRAC - 1))) : v;
            tv = tv >>> `WIN_FRAC;
            if (tv > MAXV)      rsat_win = MAXV[WL-1:0];
            else if (tv < MINV) rsat_win = MINV[WL-1:0];
            else                rsat_win = tv[WL-1:0];
        end
    endfunction

    // ---- product arithmetic ------------------------------------------------------------------------
    reg  signed [WL-1:0] f1r, f1i, f3r, f3i;
    reg  signed [WL+1:0] p12r, p12i;

    wire signed [WL-1:0] f2r = sa_q[2*WL-1:WL];
    wire signed [WL-1:0] f2i = sa_q[WL-1:0];
    wire signed [2*WL:0] mr  = f1r * f2r - f1i * f2i;           // exact real part of F1 F2
    wire signed [2*WL:0] mi  = f1r * f2i + f1i * f2r;

    function signed [WL+1:0] rsh_p;                              // round half up, shift by WL-1
        input signed [2*WL:0] v;
        reg signed [2*WL:0]   tv;
        begin
            tv = ROUND ? (v + (1 <<< (SHP - 1))) : v;
            rsh_p = tv >>> SHP;
        end
    endfunction

    wire signed [2*WL+2:0] tr = p12r * f3r + p12i * f3i;        // P12 conj(F3)
    wire signed [2*WL+2:0] ti = p12i * f3r - p12r * f3i;

    // ---- control -------------------------------------------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            st <= S_IDLE;  busy <= 1'b0;  done <= 1'b0;
            t <= 0;  k <= 0;  f1 <= 0;  f2 <= 0;  cnt <= 0;
            f_addr <= 0;  s_we <= 0;  s_waddr <= 0;  s_wdata <= 0;  ra_addr <= 0;  rb_addr <= 0;
            acc_we <= 0;  acc_waddr <= 0;  acc_ra <= 0;  accr_wdata <= 0;  acci_wdata <= 0;
            wacc_r <= 0;  wacc_i <= 0;  f1r <= 0;  f1i <= 0;  f3r <= 0;  f3i <= 0;  p12r <= 0;  p12i <= 0;
        end else begin
            done   <= 1'b0;
            s_we   <= 1'b0;
            acc_we <= 1'b0;
            case (st)
                S_IDLE: begin
                    if (clr) begin
                        st <= S_CLR;  busy <= 1'b1;  cnt <= 0;
                    end else if (go) begin
                        st <= W_A;  busy <= 1'b1;  k <= 0;  t <= 0;  wacc_r <= 0;  wacc_i <= 0;
                    end
                end
                S_CLR: begin
                    acc_we <= 1'b1;  acc_waddr <= cnt;  accr_wdata <= 0;  acci_wdata <= 0;
                    cnt <= cnt + 1'b1;
                    if (cnt == {2*LOGB{1'b1}}) st <= S_END;
                end
                // ---- window: five taps per bin -------------------------------------------------------
                W_A: begin
                    f_addr <= idxw;
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
                    s_we <= 1'b1;  s_waddr <= k;
                    s_wdata <= {rsat_win(wacc_r), rsat_win(wacc_i)};
                    wacc_r <= 0;  wacc_i <= 0;  t <= 0;
                    if (k == NBINS - 1) begin
                        f1 <= 0;
                        st <= R_A;
                    end else begin
                        k <= k + 9'd1;
                        st <= W_A;
                    end
                end
                // ---- triple loop -----------------------------------------------------------------------------
                R_A: begin
                    ra_addr <= {{(9-LOGB){1'b0}}, f1};
                    st <= R_B;
                end
                R_B: st <= R_C;
                R_C: begin
                    f1r <= sa_q[2*WL-1:WL];  f1i <= sa_q[WL-1:0];
                    f2 <= 0;
                    st <= T_A;
                end
                T_A: begin
                    ra_addr <= {{(9-LOGB){1'b0}}, f2};
                    rb_addr <= {{(9-LOGB){1'b0}}, f1} + {{(9-LOGB){1'b0}}, f2};
                    acc_ra  <= {f1, f2};
                    st <= T_B;
                end
                T_B: st <= T_C;
                T_C: begin
                    p12r <= rsh_p(mr);
                    p12i <= rsh_p(mi);
                    f3r <= sb_q[2*WL-1:WL];  f3i <= sb_q[WL-1:0];
                    st <= T_D;
                end
                T_D: begin
                    acc_we <= 1'b1;  acc_waddr <= {f1, f2};
                    accr_wdata <= accr_q + tr;
                    acci_wdata <= acci_q + ti;
                    if (f2 == MAXB - 1) begin
                        if (f1 == MAXB - 1) st <= S_END;
                        else begin f1 <= f1 + 1'b1;  st <= R_A; end
                    end else begin
                        f2 <= f2 + 1'b1;
                        st <= T_A;
                    end
                end
                default: begin                                  // S_END
                    busy <= 1'b0;  done <= 1'b1;  st <= S_IDLE;
                end
            endcase
        end
    end
endmodule
