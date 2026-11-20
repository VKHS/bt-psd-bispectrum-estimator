// Central rotational unit (CRU): y = x * exp(-j q pi/512) with a constant-coefficient rotator.
//
// Bit-exact with btpsd.cordic.rotate_fixed:
//   1. octant reduction   quad = q[9:8]; r = q[7:0]; if r > 128 { quad += 1; r = 256 - r }
//   2. ROM (a, b, inv)    for r = 0..128 (cru_rom.v, generated)
//   3. rotation           yr = xr a +/- xi b,  yi = xi a -/+ xr b   (+/-: clockwise / counter-clockwise)
//   4. gain correction    z = round(y * inv / 2^14)   (one rounding, half up, or floor if ROUND = 0)
//   5. quadrant           multiply by (-j)^quad: swap / negate, then saturate to WL bits
//
// Four register stages (LAT_CRU = 4): [octant, ROM] [rotation] [gain, shift] [quadrant, saturate].
// A new input may be presented every cycle.
`timescale 1ns / 1ps

module cru_rotate #(
    parameter WL    = 16,
    parameter ROUND = 1
) (
    input                      clk,
    input  signed [WL-1:0]     xr,
    input  signed [WL-1:0]     xi,
    input         [9:0]        q,
    output signed [WL-1:0]     yr,
    output signed [WL-1:0]     yi
);
    localparam signed [WL+12:0] MAXV = (1 <<< (WL - 1)) - 1;
    localparam signed [WL+12:0] MINV = -(1 <<< (WL - 1));
    localparam integer          RND_ADD = ROUND ? 8192 : 0;   // 2^(14-1)

    // ---- stage 0 (combinational): octant reduction and ROM ------------------------
    wire        flip   = (q[7:0] > 8'd128);
    wire [1:0]  quad_c = q[9:8] + {1'b0, flip};
    wire [7:0]  r_c    = flip ? (9'd256 - {1'b0, q[7:0]}) : q[7:0];
    wire [7:0]  a_c, b_c;
    wire [14:0] inv_c;
    cru_rom u_rom (.r(r_c), .a(a_c), .b(b_c), .inv(inv_c));

    // ---- stage 1 register ----------------------------------------------------------
    reg signed [WL-1:0] x1r, x1i;
    reg        [7:0]    a1, b1;
    reg        [14:0]   inv1;
    reg        [1:0]    quad1;
    reg                 ccw1;
    always @(posedge clk) begin
        x1r <= xr;  x1i <= xi;
        a1 <= a_c;  b1 <= b_c;  inv1 <= inv_c;
        quad1 <= quad_c;  ccw1 <= flip;
    end

    // ---- stage 2: rotation by the integer pair (a, b) ------------------------------
    wire signed [8:0]    sa = {1'b0, a1};
    wire signed [8:0]    sb = {1'b0, b1};
    wire signed [WL+8:0] p_ra = x1r * sa;
    wire signed [WL+8:0] p_ia = x1i * sa;
    wire signed [WL+8:0] p_rb = x1r * sb;
    wire signed [WL+8:0] p_ib = x1i * sb;
    wire signed [WL+9:0] yr_c = ccw1 ? (p_ra - p_ib) : (p_ra + p_ib);
    wire signed [WL+9:0] yi_c = ccw1 ? (p_ia + p_rb) : (p_ia - p_rb);

    reg signed [WL+9:0] y2r, y2i;
    reg        [14:0]   inv2;
    reg        [1:0]    quad2;
    always @(posedge clk) begin
        y2r <= yr_c;  y2i <= yi_c;
        inv2 <= inv1;  quad2 <= quad1;
    end

    // ---- stage 3: gain correction y * inv / 2^14 -----------------------------------
    wire signed [15:0]   sinv = {1'b0, inv2};
    wire signed [WL+25:0] g_r = y2r * sinv;
    wire signed [WL+25:0] g_i = y2i * sinv;
    wire signed [WL+25:0] h_r = g_r + RND_ADD;
    wire signed [WL+25:0] h_i = g_i + RND_ADD;

    reg signed [WL+11:0] z3r, z3i;
    reg        [1:0]     quad3;
    always @(posedge clk) begin
        z3r <= h_r >>> 14;
        z3i <= h_i >>> 14;
        quad3 <= quad2;
    end

    // ---- stage 4: quadrant (-j)^quad, saturation -----------------------------------
    reg signed [WL+12:0] o_r, o_i;
    always @* begin
        case (quad3)
            2'd0:    begin o_r = z3r;   o_i = z3i;   end
            2'd1:    begin o_r = z3i;   o_i = -z3r;  end
            2'd2:    begin o_r = -z3r;  o_i = -z3i;  end
            default: begin o_r = -z3i;  o_i = z3r;   end
        endcase
    end

    function signed [WL-1:0] sat;
        input signed [WL+12:0] v;
        begin
            if (v > MAXV)      sat = MAXV[WL-1:0];
            else if (v < MINV) sat = MINV[WL-1:0];
            else               sat = v[WL-1:0];
        end
    endfunction

    reg signed [WL-1:0] y4r, y4i;
    always @(posedge clk) begin
        y4r <= sat(o_r);
        y4i <= sat(o_i);
    end
    assign yr = y4r;
    assign yi = y4i;
endmodule
