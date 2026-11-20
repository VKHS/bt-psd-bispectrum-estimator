// Self-checking testbench of the memory-based 512-point FFT against btpsd.r8fft.r8fft_fixed (WL = 16, 4 chunks).
// Loads each chunk in natural order, runs the three passes, reads the spectrum in natural frequency order.
// Run from the repository root:  make -C sim fft
`timescale 1ns / 1ps
`include "fft_wl16.vh"

module tb_fft512_mem;
    localparam WL = 16, NCH = `FFT_NCHUNKS;

    reg clk = 1'b0, rst_n = 1'b0;
    always #5 clk = ~clk;

    reg [WL-1:0] in_re  [0:NCH*512-1];
    reg [WL-1:0] in_im  [0:NCH*512-1];
    reg [WL-1:0] out_re [0:NCH*512-1];
    reg [WL-1:0] out_im [0:NCH*512-1];

    reg                 ld_en, start;
    reg  [8:0]          ld_addr, rd_k;
    reg  signed [WL-1:0] ld_re, ld_im;
    wire                busy, done;
    wire signed [WL-1:0] rd_re, rd_im;

    fft512_mem #(.WL(WL), .ROUND(1)) dut (
        .clk(clk), .rst_n(rst_n), .ld_en(ld_en), .ld_addr(ld_addr), .ld_re(ld_re), .ld_im(ld_im),
        .start(start), .busy(busy), .done(done), .rd_k(rd_k), .rd_re(rd_re), .rd_im(rd_im));

    integer c, n, errs, cycles;

    // Cycle-by-cycle trace of the first transform, sampled at the rising edge like the DUT does.  Cycle 0 is the cycle in
    // which start is high; the lines (TRACE,fft,...) are analysed by python/scripts/make_sim_figures.py.
    reg     tracing = 1'b0;
    integer tcyc = 0, ttail = 0;
    always @(posedge clk) if (tracing) begin
        $display("TRACE,fft,%0d,%0d,%0d,%0d,%0d,%0d,%0d,%0d,%0d", tcyc, ld_en, start, busy, done, dut.iss_v, dut.iss_p, dut.iss_b, dut.st);
        tcyc = tcyc + 1;
        if (done) ttail = 1;
        else if (ttail != 0) ttail = ttail + 1;
        if (ttail == 6) tracing = 1'b0;
    end
    initial begin
        $readmemh("sim/vectors/fft_wl16_in_re.mem", in_re);
        $readmemh("sim/vectors/fft_wl16_in_im.mem", in_im);
        $readmemh("sim/vectors/fft_wl16_out_re.mem", out_re);
        $readmemh("sim/vectors/fft_wl16_out_im.mem", out_im);
        errs = 0;  ld_en = 1'b0;  start = 1'b0;  ld_addr = 0;  rd_k = 0;  ld_re = 0;  ld_im = 0;
        repeat (3) @(negedge clk);
        rst_n = 1'b1;
        repeat (2) @(negedge clk);
        for (c = 0; c < NCH; c = c + 1) begin
            for (n = 0; n < 512; n = n + 1) begin
                @(negedge clk);
                ld_en = 1'b1;  ld_addr = n;
                ld_re = in_re[c*512+n];  ld_im = in_im[c*512+n];
                if (c == 0 && n == 510) begin tracing = 1'b1; tcyc = -3; end
            end
            @(negedge clk);  ld_en = 1'b0;
            @(negedge clk);  start = 1'b1;
            @(negedge clk);  start = 1'b0;
            // start is sampled at the first rising edge after it is raised; the first butterfly issues in the cycle that
            // begins there, so `cycles` counts rising edges from that one to the edge at which done becomes visible.
            cycles = 0;
            while (!done && cycles < 1000) begin @(negedge clk); cycles = cycles + 1; end
            if (!done) begin
                errs = errs + 1;
                $display("chunk %0d: no done pulse within 1000 cycles", c);
            end else if (cycles != `FFT_TOTAL_CYCLES) begin
                errs = errs + 1;
                $display("chunk %0d: transform took %0d cycles, model schedule_cycles says %0d", c, cycles, `FFT_TOTAL_CYCLES);
            end else begin
                $display("chunk %0d: transform took %0d cycles from the start edge to done (= model: %0d stall cycles included)",
                         c, cycles, `FFT_STALL_CYCLES);
            end
            for (n = 0; n < 513; n = n + 1) begin
                @(negedge clk);
                if (n >= 1) begin
                    if ({rd_re, rd_im} !== {out_re[c*512+n-1], out_im[c*512+n-1]}) begin
                        errs = errs + 1;
                        if (errs <= 10) $display("MISMATCH chunk %0d bin %0d: got %h %h expected %h %h",
                                                 c, n-1, rd_re, rd_im, out_re[c*512+n-1], out_im[c*512+n-1]);
                    end
                end
                if (n < 512) rd_k = n;
            end
        end
        if (errs == 0) $display("TB_RESULT: PASS tb_fft512_mem (%0d chunks)", NCH);
        else           $display("TB_RESULT: FAIL tb_fft512_mem (%0d mismatches)", errs);
        $finish;
    end

    initial begin
        #20_000_000;
        $display("TB_RESULT: FAIL tb_fft512_mem (timeout)");
        $finish;
    end
endmodule
