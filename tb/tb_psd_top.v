// Self-checking testbench of the whole PSD estimator against btpsd.estimator.psd_fixed
// (WL = 12, Hann window, 4 chunks = 3 segments, serial sum-of-odds squarer).
// Run from the repository root:  make -C sim psd
`timescale 1ns / 1ps
`include "psd_wl12.vh"

module tb_psd_top;
    localparam WL = 12, ACCW = 2 * WL + 8, NS = `PSD_NCHUNKS * 512;

    reg clk = 1'b0, rst_n = 1'b0;
    always #5 clk = ~clk;

    reg [WL-1:0]   xmem [0:NS-1];
    reg [ACCW-1:0] pexp [0:1023];

    reg                  in_valid;
    reg  signed [WL-1:0] in_data;
    wire                 in_ready, busy;
    reg  [9:0]           res_addr;
    wire [ACCW-1:0]      res_data;

    bt_psd_top #(.WL(WL), .ROUND(1), .SQ_SERIAL(1)) dut (
        .clk(clk), .rst_n(rst_n), .avg_shift(`PSD_AVG_SHIFT),
        .in_valid(in_valid), .in_data(in_data), .in_ready(in_ready), .busy(busy),
        .res_addr(res_addr), .res_data(res_data));

    // cycles spent in the FFT phase (T_FFT = 3) and in the segment phase of the engine (T_SEG = 4) of the top level
    integer c_fft = 0, c_seg = 0;
    always @(posedge clk) if (rst_n) begin
        if (dut.st == 3'd3) c_fft = c_fft + 1;
        if (dut.st == 3'd4) c_seg = c_seg + 1;
    end

    // State occupancy (cycles per state of the top-level controller and of the segment engine, whole run) and the
    // residence time of the engine in W_S (= 4'd11, waiting for the squarer) per bin, binned by floor(log2(cycles)).
    // Printed as OCC and HIST lines, analysed by python/scripts/make_sim_figures.py.
    integer occ_top [0:7];
    integer occ_eng [0:15];
    integer ws_hist [0:15];
    integer ws_len = 0, lg, oi;
    initial begin
        for (oi = 0; oi < 8; oi = oi + 1)  occ_top[oi] = 0;
        for (oi = 0; oi < 16; oi = oi + 1) occ_eng[oi] = 0;
        for (oi = 0; oi < 16; oi = oi + 1) ws_hist[oi] = 0;
    end
    always @(posedge clk) if (rst_n) begin
        occ_top[dut.st]       = occ_top[dut.st] + 1;
        occ_eng[dut.u_eng.st] = occ_eng[dut.u_eng.st] + 1;
        if (dut.u_eng.st == 4'd11) ws_len = ws_len + 1;
        else if (ws_len != 0) begin
            lg = 0;
            while ((ws_len >> (lg + 1)) != 0) lg = lg + 1;
            ws_hist[lg] = ws_hist[lg] + 1;
            ws_len = 0;
        end
    end

    integer i, errs, cyc0, cyc;
    initial begin
        $readmemh("sim/vectors/psd_wl12_x.mem", xmem);
        $readmemh("sim/vectors/psd_wl12_P.mem", pexp);
        errs = 0;  in_valid = 1'b0;  in_data = 0;  res_addr = 0;
        repeat (3) @(negedge clk);
        rst_n = 1'b1;
        cyc0 = 0;
        i = 0;
        while (i < NS) begin
            @(negedge clk);
            if (in_ready) begin
                in_valid = 1'b1;  in_data = xmem[i];  i = i + 1;
            end else begin
                in_valid = 1'b0;
            end
        end
        @(negedge clk);  in_valid = 1'b0;
        repeat (2) @(negedge clk);
        while (!in_ready) @(negedge clk);
        $display("all %0d segments processed at t = %0d ns", `PSD_S, $time);
        $display("T_FFT %0d cycles for %0d chunks, T_SEG %0d cycles for %0d segments (serial squarer: %0d squarer cycles in the model)",
                 c_fft, `PSD_NCHUNKS, c_seg, `PSD_S, `PSD_SQUARER_CYCLES);
        for (i = 0; i < 8; i = i + 1)  if (occ_top[i] != 0) $display("OCC,psd_top,%0d,%0d", i, occ_top[i]);
        for (i = 0; i < 16; i = i + 1) if (occ_eng[i] != 0) $display("OCC,psd_engine,%0d,%0d", i, occ_eng[i]);
        for (i = 0; i < 16; i = i + 1) if (ws_hist[i] != 0) $display("HIST,psd_sq_wait_log2,%0d,%0d", i, ws_hist[i]);
        for (i = 0; i < 1025; i = i + 1) begin
            @(negedge clk);
            if (i >= 1) begin
                if (res_data !== pexp[i-1]) begin
                    errs = errs + 1;
                    if (errs <= 10) $display("MISMATCH bin %0d: got %0d expected %0d", i-1, res_data, pexp[i-1]);
                end
            end
            if (i < 1024) res_addr = i;
        end
        if (errs == 0) $display("TB_RESULT: PASS tb_psd_top (1024 bins)");
        else           $display("TB_RESULT: FAIL tb_psd_top (%0d mismatches)", errs);
        $finish;
    end

    initial begin
        #100_000_000;
        $display("TB_RESULT: FAIL tb_psd_top (timeout)");
        $finish;
    end
endmodule
