// Self-checking testbench of mag_sq_odd (serial sum-of-odds and multiplier variants), 1024 vectors.
// Run from the repository root:  make -C sim sq
`timescale 1ns / 1ps

module tb_mag_sq_odd;
    localparam WL = 16, N = 1024;

    reg clk = 1'b0, rst_n = 1'b0;
    always #5 clk = ~clk;

    reg [2*WL-1:0] vin  [0:N-1];
    reg [2*WL-1:0] vout [0:N-1];

    reg  signed [WL-1:0] re, im;
    reg                  start;
    wire                 busy_s, done_s, busy_m, done_m;
    wire [2*WL-1:0]      p_s, p_m;

    mag_sq_odd #(.WL(WL), .SERIAL(1)) dut_s (.clk(clk), .rst_n(rst_n), .start(start), .re(re), .im(im), .busy(busy_s), .done(done_s), .p(p_s));
    mag_sq_odd #(.WL(WL), .SERIAL(0)) dut_m (.clk(clk), .rst_n(rst_n), .start(start), .re(re), .im(im), .busy(busy_m), .done(done_m), .p(p_m));

    integer i, errs, cyc, maxabs, ar, ai;
    initial begin
        $readmemh("sim/vectors/sq_wl16_in.mem", vin);
        $readmemh("sim/vectors/sq_wl16_out.mem", vout);
        errs = 0;  start = 1'b0;  re = 0;  im = 0;
        repeat (3) @(negedge clk);
        rst_n = 1'b1;
        for (i = 0; i < N; i = i + 1) begin
            @(negedge clk);
            {re, im} = vin[i];
            start = 1'b1;
            @(negedge clk);
            start = 1'b0;
            cyc = 1;
            while (!done_s && cyc < 2000) begin @(negedge clk); cyc = cyc + 1; end
            ar = (re < 0) ? -re : re;
            ai = (im < 0) ? -im : im;
            maxabs = (ar > ai) ? ar : ai;
            if (p_s !== vout[i]) begin
                errs = errs + 1;
                if (errs <= 10) $display("MISMATCH serial vector %0d: got %0d expected %0d", i, p_s, vout[i]);
            end
            if (p_m !== vout[i]) begin
                errs = errs + 1;
                if (errs <= 10) $display("MISMATCH multiplier vector %0d: got %0d expected %0d", i, p_m, vout[i]);
            end
            if (cyc > maxabs + 3) begin
                errs = errs + 1;
                if (errs <= 10) $display("CYCLES vector %0d: %0d cycles for max(|re|,|im|) = %0d", i, cyc, maxabs);
            end
        end
        if (errs == 0) $display("TB_RESULT: PASS tb_mag_sq_odd (%0d vectors)", N);
        else           $display("TB_RESULT: FAIL tb_mag_sq_odd (%0d errors)", errs);
        $finish;
    end
endmodule
