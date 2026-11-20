// Self-checking testbench of the direct bispectrum estimator against btpsd.bispectrum.bispectrum_fixed
// (WL = 12, Hann window, 5 frames with hop 256, 32 x 32 bins).
// Run from the repository root:  make -C sim bsp
`timescale 1ns / 1ps
`include "bsp_wl12.vh"

module tb_bispectrum_top;
    localparam WL = 12, LOGB = 5, MAXB = 1 << LOGB, BACCW = 2 * WL + 10, NS = `BSP_NSAMPLES;

    reg clk = 1'b0, rst_n = 1'b0;
    always #5 clk = ~clk;

    reg [WL-1:0]    xmem [0:NS-1];
    reg [BACCW-1:0] bre  [0:MAXB*MAXB-1];
    reg [BACCW-1:0] bim  [0:MAXB*MAXB-1];

    reg                   in_valid;
    reg  signed [WL-1:0]  in_data;
    wire                  in_ready, busy;
    reg  [2*LOGB-1:0]     res_idx;
    wire signed [BACCW-1:0] res_re, res_im;

    bt_bispectrum_top #(.WL(WL), .ROUND(1), .LOGB(LOGB)) dut (
        .clk(clk), .rst_n(rst_n), .avg_shift(`BSP_AVG_SHIFT),
        .in_valid(in_valid), .in_data(in_data), .in_ready(in_ready), .busy(busy),
        .res_idx(res_idx), .res_re(res_re), .res_im(res_im));

    integer i, errs;
    initial begin
        if (`BSP_MAX_BIN != MAXB) begin
            $display("TB_RESULT: FAIL tb_bispectrum_top (vectors are for MAX_BIN = %0d, testbench for %0d)", `BSP_MAX_BIN, MAXB);
            $finish;
        end
        $readmemh("sim/vectors/bsp_wl12_x.mem", xmem);
        $readmemh("sim/vectors/bsp_wl12_B_re.mem", bre);
        $readmemh("sim/vectors/bsp_wl12_B_im.mem", bim);
        errs = 0;  in_valid = 1'b0;  in_data = 0;  res_idx = 0;
        repeat (3) @(negedge clk);
        rst_n = 1'b1;
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
        $display("all %0d frames processed at t = %0d ns", `BSP_S, $time);
        for (i = 0; i < MAXB*MAXB + 1; i = i + 1) begin
            @(negedge clk);
            if (i >= 1) begin
                if ({res_re, res_im} !== {bre[i-1], bim[i-1]}) begin
                    errs = errs + 1;
                    if (errs <= 10) $display("MISMATCH f1 = %0d f2 = %0d: got %0d %0d expected %0d %0d",
                                             (i-1) / MAXB, (i-1) % MAXB, res_re, res_im, $signed(bre[i-1]), $signed(bim[i-1]));
                end
            end
            if (i < MAXB*MAXB) res_idx = i;
        end
        if (errs == 0) $display("TB_RESULT: PASS tb_bispectrum_top (%0d bins)", MAXB*MAXB);
        else           $display("TB_RESULT: FAIL tb_bispectrum_top (%0d mismatches)", errs);
        $finish;
    end

    initial begin
        #50_000_000;
        $display("TB_RESULT: FAIL tb_bispectrum_top (timeout)");
        $finish;
    end
endmodule
