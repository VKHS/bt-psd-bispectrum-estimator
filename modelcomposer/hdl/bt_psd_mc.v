// Model Composer / System Generator Black Box wrapper of the chunk-merge PSD estimator (bt_psd_top).
//
// Why a wrapper.  bt_psd_top takes its samples through a ready/valid handshake and holds the producer off for the FFT and the
// segment merge (hundreds of cycles per chunk).  A Simulink signal source cannot be held off: it emits one sample per sample
// period.  The wrapper therefore RECORDS the first NCHUNKS * 512 samples that arrive (in_valid high) in a RAM, REPLAYS them into
// bt_psd_top at the pace the core accepts them, and raises done when the accumulated power spectrum is complete.
//
//   rst --> CAPTURE --(NSAMP samples seen)--> PRIME --> REPLAY --(last sample accepted)--> DRAIN --(core idle)--> DONE
//            in_valid & ce write the RAM       1 cycle   core in_valid, core in_ready       core busy              rd_addr/rd_data
//
// A one-cycle pulse on rst (synchronous, active high) re-arms the block; the estimator core is reset together with it and clears
// its accumulators by itself (about 1 024 cycles, hidden behind the capture phase).
//
// Ports
//   clk, ce       clock and clock enable.  Single-rate use: ce is expected to be 1 in every cycle.  It qualifies in_valid only; the
//                 replay and the estimator run on every clk edge, so the block must not be used in a multi-rate design.
//   rst           synchronous reset, active high
//   in_valid      1 for every sample of the record; samples that arrive while the block is not capturing are ignored
//   in_data       sample, WL-bit two's complement (Fix_WL_(WL-1) in Simulink: -1 <= x < 1)
//   rd_addr       bin to read, 0 .. 1023 (two-sided spectrum); rd_data follows one clock later
//   busy          1 from the first captured sample until done
//   done          1 while the spectrum is complete and readable (stays until the next rst)
//   rd_data       accumulated, shifted power of bin rd_addr (ACCW = 2 WL + 8 bit, unsigned); valid while done = 1.  Scale to the
//                 Welch estimate with btpsd.estimator.hw_power_scale(WL, NCHUNKS - 1, avg_shift).
//
// The number of segments is NCHUNKS - 1 and avg_shift = ceil(log2(NCHUNKS - 1)) is computed here, so the power spectrum is the
// average over the segments exactly when NCHUNKS - 1 is a power of two and scaled by (NCHUNKS - 1) / 2^avg_shift otherwise.
//
// Bit-exact with btpsd.estimator.psd_fixed on the first NCHUNKS * 512 samples (tb/tb_bt_psd_mc.v).
`timescale 1ns / 1ps
`include "btpsd_consts.vh"

module bt_psd_mc #(
    parameter WL        = 12,
    parameter ROUND     = 1,
    parameter SQ_SERIAL = 1,
    parameter NCHUNKS   = 5,
    parameter ACCW      = 2 * WL + 8,
    parameter signed [15:0] WM2 = `WIN_HANN_M2,
    parameter signed [15:0] WM1 = `WIN_HANN_M1,
    parameter signed [15:0] W0  = `WIN_HANN_0,
    parameter signed [15:0] WP1 = `WIN_HANN_P1,
    parameter signed [15:0] WP2 = `WIN_HANN_P2
) (
    input                      clk,
    input                      ce,
    input                      rst,
    input                      in_valid,
    input  signed [WL-1:0]     in_data,
    input      [9:0]           rd_addr,
    output                     busy,
    output                     done,
    output     [ACCW-1:0]      rd_data
);
    function integer clog2;
        input integer v;
        integer i;
        begin
            clog2 = 0;
            for (i = v - 1; i > 0; i = i >> 1) clog2 = clog2 + 1;
        end
    endfunction

    localparam integer NSAMP = NCHUNKS * 512;
    localparam integer AW    = clog2(NSAMP);
    localparam [4:0]   AVG_SHIFT = clog2(NCHUNKS - 1);

    localparam [2:0] S_CAPTURE = 3'd0, S_PRIME = 3'd1, S_REPLAY = 3'd2, S_DRAIN = 3'd3, S_DONE = 3'd4;

    reg [2:0]    state  = S_CAPTURE;
    reg [AW:0]   wr_cnt = 0;                 // samples captured
    reg [AW:0]   rp     = 0;                 // replay pointer: index of the sample presented to the core
    reg [WL-1:0] ram [0:NSAMP-1];
    reg [WL-1:0] rdata;

    // ---- reset of the estimator core: low from one cycle after rst until three cycles after rst falls ----------------------
    reg [3:0] rst_sr = 4'b0000;
    always @(posedge clk) rst_sr <= {rst_sr[2:0], ~rst};
    wire core_rst_n = &rst_sr;

    // ---- capture ---------------------------------------------------------------------------------------------------------------
    wire cap = (state == S_CAPTURE) && ce && in_valid;
    always @(posedge clk) if (cap) ram[wr_cnt[AW-1:0]] <= in_data;

    // ---- replay: synchronous read, the word presented in a cycle is ram[rp] ----------------------------------------------------
    wire         core_in_valid = (state == S_REPLAY);
    wire         core_in_ready, core_busy;
    wire         accept  = core_in_valid && core_in_ready;
    wire [AW:0]  rp_next = accept ? rp + 1'b1 : rp;
    wire [AW-1:0] ra     = (rp_next >= NSAMP) ? {AW{1'b0}} : rp_next[AW-1:0];
    always @(posedge clk) rdata <= ram[ra];

    always @(posedge clk) begin
        if (rst) begin
            state <= S_CAPTURE;  wr_cnt <= 0;  rp <= 0;
        end else begin
            rp <= rp_next;
            case (state)
                S_CAPTURE: if (cap) begin
                    wr_cnt <= wr_cnt + 1'b1;
                    if (wr_cnt == NSAMP - 1) state <= S_PRIME;
                end
                S_PRIME:   state <= S_REPLAY;
                S_REPLAY:  if (accept && rp == NSAMP - 1) state <= S_DRAIN;
                S_DRAIN:   if (!core_busy) state <= S_DONE;      // busy is already 1 in the first cycle after the last sample
                default:   state <= S_DONE;
            endcase
        end
    end

    assign busy = (state == S_CAPTURE) ? (wr_cnt != 0) : (state != S_DONE);
    assign done = (state == S_DONE);

    bt_psd_top #(
        .WL(WL), .ROUND(ROUND), .SQ_SERIAL(SQ_SERIAL), .ACCW(ACCW),
        .WM2(WM2), .WM1(WM1), .W0(W0), .WP1(WP1), .WP2(WP2)
    ) u_core (
        .clk(clk), .rst_n(core_rst_n), .avg_shift(AVG_SHIFT),
        .in_valid(core_in_valid), .in_data(rdata), .in_ready(core_in_ready), .busy(core_busy),
        .res_addr(rd_addr), .res_data(rd_data)
    );
endmodule
