// Self-checking testbench of r8_butterfly against dft8_fixed + rotate_fixed (WL = 16, 2048 vectors).
// Run from the repository root:  make -C sim bfly
`timescale 1ns / 1ps
`include "btpsd_params.vh"

module tb_r8_butterfly;
    localparam WL = 16, N = 2048, LAT = `LAT_BFLY;

    reg clk = 1'b0;
    always #5 clk = ~clk;

    reg [16*WL+69:0] vin  [0:N-1];
    reg [16*WL-1:0]  vout [0:N-1];

    reg  [8*WL-1:0] xr, xi;
    reg  [69:0]     tq;
    wire [8*WL-1:0] yr, yi;

    r8_butterfly #(.WL(WL), .ROUND(1)) dut (.clk(clk), .xr(xr), .xi(xi), .tq(tq), .yr(yr), .yi(yi));

    integer i, errs;
    initial begin
        $readmemh("sim/vectors/bfly_wl16_in.mem", vin);
        $readmemh("sim/vectors/bfly_wl16_out.mem", vout);
        errs = 0;
        xr = 0;  xi = 0;  tq = 0;
        for (i = 0; i < N + LAT; i = i + 1) begin
            @(negedge clk);
            if (i >= LAT) begin
                if ({yr, yi} !== vout[i-LAT]) begin
                    errs = errs + 1;
                    if (errs <= 5) $display("MISMATCH vector %0d:\n  got      %h\n  expected %h", i-LAT, {yr, yi}, vout[i-LAT]);
                end
            end
            if (i < N) {xr, xi, tq} = vin[i];
        end
        if (errs == 0) $display("TB_RESULT: PASS tb_r8_butterfly (%0d vectors)", N);
        else           $display("TB_RESULT: FAIL tb_r8_butterfly (%0d mismatches)", errs);
        $finish;
    end
endmodule
