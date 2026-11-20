// Radix-8 butterfly: scaled (1/8) 8-point DFT followed by the output twiddle rotations.
//
// Bit-exact with  btpsd.r8fft.dft8_fixed  followed by  btpsd.cordic.rotate_fixed  (output k rotated by tq[k]):
//   stage 1  distance 4: s_i = (x_i + x_{i+4}) >> 1,  d_i = (x_i - x_{i+4}) >> 1,  d_i *= W8^i (CRU, q = 0, 128, 256, 384)
//   stage 2  distance 2: c0 = g0 + g2, c2 = g0 - g2, c1 = g1 + g3, c3 = (-j)(g1 - g3)          (each result >> 1)
//   stage 3  distance 1: out[slot] = c0 + c1, out[slot+4] = c0 - c1, out[slot+2] = c2 + c3, ... (each result >> 1)
//   twiddle  Y[k] *= exp(-j tq[k] pi/512), k = 1..7 (CRU); Y[0] is passed through.
// Every add/subtract is followed by one rounded right shift and clamped to WL bits (only the rounded difference of the
// two extreme codes can exceed the range, by one LSB); the (-j) negation saturates.
//
// Ports are packed buses: slot m of x/y at bits [m*WL +: WL]; twiddle slot k = 1..7 at tq[(k-1)*10 +: 10].
// tq is presented together with x.  Latency LAT_BFLY = 1 + LAT_CRU + 1 + LAT_CRU = 10 cycles; one butterfly per cycle.
`timescale 1ns / 1ps
`include "btpsd_params.vh"

module r8_butterfly #(
    parameter WL    = 16,
    parameter ROUND = 1
) (
    input                  clk,
    input  [8*WL-1:0]      xr,
    input  [8*WL-1:0]      xi,
    input  [69:0]          tq,
    output [8*WL-1:0]      yr,
    output [8*WL-1:0]      yi
);
    localparam LC = `LAT_CRU;

    // ---- arithmetic helpers --------------------------------------------------------
    // Clamp a (WL+2)-bit two's-complement value to WL bits: it fits exactly when its top three bits agree.
    function signed [WL-1:0] sat_wl;
        input signed [WL+1:0] t;
        begin
            if (t[WL+1] == t[WL] && t[WL] == t[WL-1]) sat_wl = t[WL-1:0];
            else if (t[WL+1])                         sat_wl = {1'b1, {(WL-1){1'b0}}};
            else                                      sat_wl = {1'b0, {(WL-1){1'b1}}};
        end
    endfunction

    // (a + b) >> 1 and (a - b) >> 1 with rounding, clamped to WL bits.  The clamp only matters for the difference
    // of the two extreme codes: (2^(WL-1) - 1) - (-2^(WL-1)) + 1 halves to +2^(WL-1), one above the largest code.
    function signed [WL-1:0] hadd;
        input signed [WL-1:0] a;
        input signed [WL-1:0] b;
        reg signed [WL+1:0] t;
        begin
            t = a + b;
            if (ROUND) t = t + 1;
            t = t >>> 1;
            hadd = sat_wl(t);
        end
    endfunction

    function signed [WL-1:0] hsub;
        input signed [WL-1:0] a;
        input signed [WL-1:0] b;
        reg signed [WL+1:0] t;
        begin
            t = a - b;
            if (ROUND) t = t + 1;
            t = t >>> 1;
            hsub = sat_wl(t);
        end
    endfunction

    // -a, saturating the single overflow case a = -2^(WL-1)
    function signed [WL-1:0] nsat;
        input signed [WL-1:0] a;
        begin
            if (a == {1'b1, {(WL-1){1'b0}}}) nsat = {1'b0, {(WL-1){1'b1}}};
            else                             nsat = -a;
        end
    endfunction

    integer i;

    // ---- stage 1: distance-4 add/subtract -------------------------------------------
    reg [4*WL-1:0] s1r, s1i, d1r, d1i;
    always @(posedge clk) begin
        for (i = 0; i < 4; i = i + 1) begin
            s1r[i*WL +: WL] <= hadd(xr[i*WL +: WL], xr[(i+4)*WL +: WL]);
            s1i[i*WL +: WL] <= hadd(xi[i*WL +: WL], xi[(i+4)*WL +: WL]);
            d1r[i*WL +: WL] <= hsub(xr[i*WL +: WL], xr[(i+4)*WL +: WL]);
            d1i[i*WL +: WL] <= hsub(xi[i*WL +: WL], xi[(i+4)*WL +: WL]);
        end
    end

    // ---- W8^i on the difference branch (i = 1, 2, 3 through the CRU; i = 0 is delayed) -
    wire [WL-1:0] w1r, w1i, w2r, w2i, w3r, w3i;
    cru_rotate #(.WL(WL), .ROUND(ROUND)) u_w1 (.clk(clk), .xr(d1r[1*WL +: WL]), .xi(d1i[1*WL +: WL]), .q(10'd128), .yr(w1r), .yi(w1i));
    cru_rotate #(.WL(WL), .ROUND(ROUND)) u_w2 (.clk(clk), .xr(d1r[2*WL +: WL]), .xi(d1i[2*WL +: WL]), .q(10'd256), .yr(w2r), .yi(w2i));
    cru_rotate #(.WL(WL), .ROUND(ROUND)) u_w3 (.clk(clk), .xr(d1r[3*WL +: WL]), .xi(d1i[3*WL +: WL]), .q(10'd384), .yr(w3r), .yi(w3i));

    wire [4*WL-1:0] s2r, s2i;
    wire [WL-1:0]   d20r, d20i;
    delay_line #(.W(10*WL), .N(LC)) u_dl_s (
        .clk(clk),
        .d({s1r, s1i, d1r[WL-1:0], d1i[WL-1:0]}),
        .q({s2r, s2i, d20r, d20i}));

    wire [4*WL-1:0] gdr = {w3r, w2r, w1r, d20r};   // rotated difference branch, slot m at m*WL
    wire [4*WL-1:0] gdi = {w3i, w2i, w1i, d20i};

    // ---- stages 2 and 3: distance-2 and distance-1 add/subtract (one register stage) ---
    reg [8*WL-1:0]  y3r, y3i;
    reg [4*WL-1:0]  gr, gi;
    reg signed [WL-1:0] c0r, c0i, c1r, c1i, c2r, c2i, c3r, c3i, e3r, e3i;
    integer grp;
    always @(posedge clk) begin
        for (grp = 0; grp < 2; grp = grp + 1) begin
            gr = grp ? gdr : s2r;
            gi = grp ? gdi : s2i;
            c0r = hadd(gr[0*WL +: WL], gr[2*WL +: WL]);
            c0i = hadd(gi[0*WL +: WL], gi[2*WL +: WL]);
            c2r = hsub(gr[0*WL +: WL], gr[2*WL +: WL]);
            c2i = hsub(gi[0*WL +: WL], gi[2*WL +: WL]);
            c1r = hadd(gr[1*WL +: WL], gr[3*WL +: WL]);
            c1i = hadd(gi[1*WL +: WL], gi[3*WL +: WL]);
            e3r = hsub(gr[1*WL +: WL], gr[3*WL +: WL]);
            e3i = hsub(gi[1*WL +: WL], gi[3*WL +: WL]);
            c3r = e3i;            // (-j)(e3r + j e3i) = e3i - j e3r
            c3i = nsat(e3r);
            y3r[(0+grp)*WL +: WL] <= hadd(c0r, c1r);
            y3i[(0+grp)*WL +: WL] <= hadd(c0i, c1i);
            y3r[(4+grp)*WL +: WL] <= hsub(c0r, c1r);
            y3i[(4+grp)*WL +: WL] <= hsub(c0i, c1i);
            y3r[(2+grp)*WL +: WL] <= hadd(c2r, c3r);
            y3i[(2+grp)*WL +: WL] <= hadd(c2i, c3i);
            y3r[(6+grp)*WL +: WL] <= hsub(c2r, c3r);
            y3i[(6+grp)*WL +: WL] <= hsub(c2i, c3i);
        end
    end

    // ---- output twiddles --------------------------------------------------------------
    wire [69:0] tqd;
    delay_line #(.W(70), .N(LC + 2)) u_dl_tq (.clk(clk), .d(tq), .q(tqd));

    delay_line #(.W(2*WL), .N(LC)) u_dl_y0 (
        .clk(clk),
        .d({y3r[WL-1:0], y3i[WL-1:0]}),
        .q({yr[WL-1:0], yi[WL-1:0]}));

    genvar k;
    generate
        for (k = 1; k < 8; k = k + 1) begin : g_tw
            cru_rotate #(.WL(WL), .ROUND(ROUND)) u_tw (
                .clk(clk),
                .xr(y3r[k*WL +: WL]), .xi(y3i[k*WL +: WL]),
                .q(tqd[(k-1)*10 +: 10]),
                .yr(yr[k*WL +: WL]),  .yi(yi[k*WL +: WL]));
        end
    endgenerate
endmodule
