// Self-checking testbench of the in-fabric self-test bt_psd_selftest (sig_gen -> bt_psd_mc -> bt_psd_top) against the golden spectrum
// of the first 4 chunks of the generator (sim/vectors/selftest_wl12_P.mem).  One run after a reset pulse; the result is read after done.
// Run from the repository root:  make -C sim st
`timescale 1ns / 1ps
`include "selftest_wl12.vh"

module tb_bt_psd_selftest;
    localparam WL = 12, ACCW = 2 * WL + 8;

    reg clk = 1'b0;
    always #5 clk = ~clk;

    reg [ACCW-1:0] pexp [0:1023];

    reg                rst = 1'b0, noise_en = 1'b1;
    reg  [31:0]        pinc = `ST_PINC;
    reg  [3:0]         noise_shift = `ST_NOISE_SHIFT;
    reg  [9:0]         rd_addr = 0;
    wire               busy, done;
    wire [ACCW-1:0]    rd_data;

    bt_psd_selftest #(.WL(WL), .NCHUNKS(`ST_NCHUNKS)) dut (
        .clk(clk), .rst(rst), .pinc(pinc), .noise_en(noise_en), .noise_shift(noise_shift), .rd_addr(rd_addr),
        .busy(busy), .done(done), .rd_data(rd_data));

    integer i, errs, cyc;
    initial begin
        $readmemh("sim/vectors/selftest_wl12_P.mem", pexp);
        errs = 0;  cyc = 0;
        repeat (6) @(negedge clk);
        @(negedge clk);  rst = 1'b1;
        @(negedge clk);  rst = 1'b0;
        while (!done && cyc < 2_000_000) begin
            @(negedge clk);  cyc = cyc + 1;
        end
        if (!done) begin
            errs = errs + 1;  $display("done not reached");
        end
        $display("MC,selftest,cycles_from_reset_to_done,%0d", cyc);
        for (i = 0; i < 1025; i = i + 1) begin
            @(negedge clk);
            if (i >= 1 && rd_data !== pexp[i-1]) begin
                errs = errs + 1;
                if (errs <= 10) $display("MISMATCH bin %0d: got %0d expected %0d", i-1, rd_data, pexp[i-1]);
            end
            if (i < 1024) rd_addr = i;
        end
        if (errs == 0) $display("TB_RESULT: PASS tb_bt_psd_selftest (1024 bins)");
        else           $display("TB_RESULT: FAIL tb_bt_psd_selftest (%0d mismatches)", errs);
        $finish;
    end

    initial begin
        #100_000_000;
        $display("TB_RESULT: FAIL tb_bt_psd_selftest (timeout)");
        $finish;
    end
endmodule
