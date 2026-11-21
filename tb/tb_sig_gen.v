// Self-checking testbench of the test-signal generator sig_gen against btpsd.siggen.generate (sim/vectors/siggen_wl12_x.mem).
// Four runs, each after a reset:  A tone + noise (noise_shift = SG_SHIFT_A), en = 1 in every cycle;
//                                 B the same with en low in every third cycle (all registers must hold),
//                                 C noise at full scale (noise_shift = SG_SHIFT_C): the saturation of the sum,
//                                 D noise disabled (tone only).
// valid must be low until two enabled clock edges have passed and high afterwards.
// Run from the repository root:  make -C sim sg
`timescale 1ns / 1ps
`include "siggen_wl12.vh"

module tb_sig_gen;
    localparam WL = 12, N = `SG_N;

    reg clk = 1'b0;
    always #5 clk = ~clk;

    reg [WL-1:0] xexp [0:3*N-1];

    reg                 rst_n = 1'b0, en = 1'b0, noise_en = 1'b1;
    reg  [31:0]         pinc = `SG_PINC;
    reg  [3:0]          noise_shift = `SG_SHIFT_A;
    wire signed [WL-1:0] x;
    wire                valid;

    sig_gen #(.WL(WL)) dut (
        .clk(clk), .rst_n(rst_n), .en(en), .pinc(pinc), .noise_en(noise_en), .noise_shift(noise_shift),
        .x(x), .valid(valid));

    integer errs, got, ne, cyc, gapped, base;
    reg en_used;

    task run_one;
        begin
            @(negedge clk);  rst_n = 1'b0;  en = 1'b0;
            @(negedge clk);
            @(negedge clk);  rst_n = 1'b1;
            got = 0;  ne = 0;  cyc = 0;
            while (got < N && cyc < 100000) begin
                en_used = gapped ? ((cyc % 3) != 1) : 1'b1;
                en = en_used;
                @(negedge clk);                      // the clock edge in between has seen en
                cyc = cyc + 1;
                if (en_used) ne = ne + 1;
                if (valid !== (ne >= 2)) begin
                    errs = errs + 1;
                    if (errs <= 10) $display("valid = %b after %0d enabled edges", valid, ne);
                end
                if (en_used && valid === 1'b1) begin
                    if (x !== xexp[base + got]) begin
                        errs = errs + 1;
                        if (errs <= 10) $display("MISMATCH run base %0d sample %0d: got %0d expected %0d", base, got, x, $signed(xexp[base + got]));
                    end
                    got = got + 1;
                end
            end
            if (got != N) begin
                errs = errs + 1;  $display("run base %0d: only %0d of %0d samples", base, got, N);
            end
            en = 1'b0;
        end
    endtask

    initial begin
        $readmemh("sim/vectors/siggen_wl12_x.mem", xexp);
        errs = 0;
        // run A
        gapped = 0;  base = 0;  noise_en = 1'b1;  noise_shift = `SG_SHIFT_A;  run_one;
        // run B
        gapped = 1;  base = 0;  run_one;
        // run C
        gapped = 0;  base = N;  noise_shift = `SG_SHIFT_C;  run_one;
        // run D
        gapped = 0;  base = 2 * N;  noise_en = 1'b0;  noise_shift = `SG_SHIFT_A;  run_one;
        if (errs == 0) $display("TB_RESULT: PASS tb_sig_gen (4 runs x %0d samples)", N);
        else           $display("TB_RESULT: FAIL tb_sig_gen (%0d mismatches)", errs);
        $finish;
    end

    initial begin
        #20_000_000;
        $display("TB_RESULT: FAIL tb_sig_gen (timeout)");
        $finish;
    end
endmodule
