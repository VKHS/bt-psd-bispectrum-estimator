// Self-checking testbench of the Model Composer wrapper bt_psd_mc against the golden PSD vectors (WL = 12, Hann, 4 chunks = 3 segments;
// the same record and result as tb_psd_top).  Two runs of the record:
//   run 0  one sample in every cycle; junk samples are driven while the result is read out (they must be ignored);
//   run 1  after a one-cycle rst pulse, a sample only every third cycle, junk (in_valid low) in between.
// Checks: done is low after rst and high at the end, busy is seen, all 1 024 bins equal the model, re-arming works.
// Run from the repository root:  make -C sim mc      (wrapper and RTL sources)
//                                make -C sim mcb     (the include-free bundle modelcomposer/hdl/bt_psd_mc_bundle.v alone)
`timescale 1ns / 1ps
`include "psd_wl12.vh"

module tb_bt_psd_mc;
    localparam WL = 12, ACCW = 2 * WL + 8, NCH = `PSD_NCHUNKS, NS = NCH * 512;

    reg clk = 1'b0;
    always #5 clk = ~clk;

    reg [WL-1:0]   xmem [0:NS-1];
    reg [ACCW-1:0] pexp [0:1023];

    reg                  ce = 1'b1, rst = 1'b0, in_valid = 1'b0;
    reg  signed [WL-1:0] in_data = 0;
    reg  [9:0]           rd_addr = 0;
    wire                 busy, done;
    wire [ACCW-1:0]      rd_data;

    bt_psd_mc #(.WL(WL), .NCHUNKS(NCH)) dut (
        .clk(clk), .ce(ce), .rst(rst), .in_valid(in_valid), .in_data(in_data), .rd_addr(rd_addr),
        .busy(busy), .done(done), .rd_data(rd_data));

    integer i, run, errs, cyc, saw_busy, busy_before_done;
    initial begin
        $readmemh("sim/vectors/psd_wl12_x.mem", xmem);
        $readmemh("sim/vectors/psd_wl12_P.mem", pexp);
        errs = 0;
        repeat (6) @(negedge clk);
        for (run = 0; run < 2; run = run + 1) begin
            @(negedge clk);  rst = 1'b1;
            @(negedge clk);  rst = 1'b0;
            @(negedge clk);
            if (done !== 1'b0 || busy !== 1'b0) begin
                errs = errs + 1;  $display("run %0d: done=%b busy=%b after rst (expected 0, 0)", run, done, busy);
            end
            saw_busy = 0;
            cyc = 0;
            for (i = 0; i < NS; i = i + 1) begin
                if (run == 1) begin
                    in_valid = 1'b0;  in_data = ~xmem[i];  @(negedge clk);  cyc = cyc + 1;     // two idle cycles, junk on the bus
                    in_valid = 1'b0;  in_data = ~xmem[i];  @(negedge clk);  cyc = cyc + 1;
                end
                in_valid = 1'b1;  in_data = xmem[i];  @(negedge clk);  cyc = cyc + 1;
                if (busy) saw_busy = 1;
            end
            in_valid = 1'b0;
            busy_before_done = 0;
            while (!done && cyc < 2_000_000) begin
                @(negedge clk);  cyc = cyc + 1;
                if (busy && !done) busy_before_done = 1;
            end
            if (!done) begin
                errs = errs + 1;  $display("run %0d: done not reached", run);
            end
            if (!saw_busy || !busy_before_done) begin
                errs = errs + 1;  $display("run %0d: busy not seen (capture %0d, processing %0d)", run, saw_busy, busy_before_done);
            end
            $display("MC,run,%0d,cycles_from_start_to_done,%0d", run, cyc);
            for (i = 0; i < 1025; i = i + 1) begin
                @(negedge clk);
                if (run == 0) begin
                    in_valid = (i % 5) != 0;  in_data = i;                                     // junk after done: must be ignored
                end else begin
                    in_valid = 1'b0;
                end
                if (!done) begin
                    errs = errs + 1;  if (errs <= 10) $display("run %0d: done dropped during readout (bin %0d)", run, i);
                end
                if (i >= 1) begin
                    if (rd_data !== pexp[i-1]) begin
                        errs = errs + 1;
                        if (errs <= 10) $display("MISMATCH run %0d bin %0d: got %0d expected %0d", run, i-1, rd_data, pexp[i-1]);
                    end
                end
                if (i < 1024) rd_addr = i;
            end
            in_valid = 1'b0;
        end
        if (errs == 0) $display("TB_RESULT: PASS tb_bt_psd_mc (2 runs x 1024 bins)");
        else           $display("TB_RESULT: FAIL tb_bt_psd_mc (%0d mismatches)", errs);
        $finish;
    end

    initial begin
        #200_000_000;
        $display("TB_RESULT: FAIL tb_bt_psd_mc (timeout)");
        $finish;
    end
endmodule
