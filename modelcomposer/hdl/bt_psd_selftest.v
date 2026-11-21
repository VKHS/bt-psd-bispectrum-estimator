// In-fabric self-test of the estimator: the test-signal generator sig_gen (tone + noise) feeds the record-and-replay wrapper
// bt_psd_mc, so the whole chain runs without Simulink and without input data.  A one-cycle pulse on rst (synchronous, active high)
// restarts generator and estimator; the estimator also starts by itself at power-up.  pinc, noise_en and noise_shift must stay
// constant during a run.  When done = 1, rd_addr -> rd_data (one clock later) reads the 1024-bin spectrum of the first
// NCHUNKS * 512 generator samples; btpsd.siggen / btpsd.estimator.psd_fixed give the expected words (tb/tb_bt_psd_selftest.v).
`timescale 1ns / 1ps

module bt_psd_selftest #(
    parameter WL        = 12,
    parameter ROUND     = 1,
    parameter SQ_SERIAL = 1,
    parameter NCHUNKS   = 5,
    parameter ACCW      = 2 * WL + 8
) (
    input                  clk,
    input                  rst,
    input      [31:0]      pinc,
    input                  noise_en,
    input      [3:0]       noise_shift,
    input      [9:0]       rd_addr,
    output                 busy,
    output                 done,
    output     [ACCW-1:0]  rd_data
);
    wire signed [WL-1:0] x;
    wire                 valid;

    sig_gen #(.WL(WL)) u_gen (
        .clk(clk), .rst_n(~rst), .en(1'b1), .pinc(pinc), .noise_en(noise_en), .noise_shift(noise_shift),
        .x(x), .valid(valid));

    bt_psd_mc #(.WL(WL), .ROUND(ROUND), .SQ_SERIAL(SQ_SERIAL), .NCHUNKS(NCHUNKS), .ACCW(ACCW)) u_mc (
        .clk(clk), .ce(1'b1), .rst(rst), .in_valid(valid), .in_data(x), .rd_addr(rd_addr),
        .busy(busy), .done(done), .rd_data(rd_data));
endmodule
