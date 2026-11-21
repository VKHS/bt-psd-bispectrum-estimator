// Self-checking testbench of the direct bispectrum estimator against btpsd.bispectrum.bispectrum_fixed
// (default: WL = 12, Hann window, 5 frames with hop 256, 32 x 32 bins).
// Run from the repository root:  make -C sim bsp
//
// Other configurations are run by python/scripts/run_rtl_matrix.py (see tb_psd_top.v): TB_VH names the generated include
// file with the BSP_* counts, which may define TB_WL, TB_PARAMS, TB_XFILE, TB_REFILE and TB_IMFILE.
`timescale 1ns / 1ps
`ifndef TB_VH
`define TB_VH "bsp_wl12.vh"
`endif
`include `TB_VH
`ifndef TB_WL
`define TB_WL 12
`endif
`ifndef TB_PARAMS
`define TB_PARAMS .WL(`TB_WL), .ROUND(1), .LOGB(5)
`endif
`ifndef TB_XFILE
`define TB_XFILE "sim/vectors/bsp_wl12_x.mem"
`endif
`ifndef TB_REFILE
`define TB_REFILE "sim/vectors/bsp_wl12_B_re.mem"
`endif
`ifndef TB_IMFILE
`define TB_IMFILE "sim/vectors/bsp_wl12_B_im.mem"
`endif

module tb_bispectrum_top;
    localparam WL = `TB_WL, LOGB = 5, MAXB = 1 << LOGB, BACCW = 2 * WL + 10, NS = `BSP_NSAMPLES;

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

    bt_bispectrum_top #(`TB_PARAMS) dut (
        .clk(clk), .rst_n(rst_n), .avg_shift(`BSP_AVG_SHIFT),
        .in_valid(in_valid), .in_data(in_data), .in_ready(in_ready), .busy(busy),
        .res_idx(res_idx), .res_re(res_re), .res_im(res_im));

    // cycles spent in the FFT phase (B_FFT = 4) and in the bispectrum engine (B_ENG = 5) of the top level
    integer c_fft = 0, c_eng = 0;
    always @(posedge clk) if (rst_n) begin
        if (dut.st == 3'd4) c_fft = c_fft + 1;
        if (dut.st == 3'd5) c_eng = c_eng + 1;
    end

    // State occupancy (cycles per state of the top-level controller and of the bispectrum engine, whole run); printed as
    // OCC lines, analysed by python/scripts/make_sim_figures.py.
    integer occ_top [0:7];
    integer occ_eng [0:15];
    integer oi;
    initial begin
        for (oi = 0; oi < 8; oi = oi + 1)  occ_top[oi] = 0;
        for (oi = 0; oi < 16; oi = oi + 1) occ_eng[oi] = 0;
    end
    always @(posedge clk) if (rst_n) begin
        occ_top[dut.st]       = occ_top[dut.st] + 1;
        occ_eng[dut.u_eng.st] = occ_eng[dut.u_eng.st] + 1;
    end

    integer i, errs;
    initial begin
        if (`BSP_MAX_BIN != MAXB) begin
            $display("TB_RESULT: FAIL tb_bispectrum_top (vectors are for MAX_BIN = %0d, testbench for %0d)", `BSP_MAX_BIN, MAXB);
            $finish;
        end
        $readmemh(`TB_XFILE, xmem);
        $readmemh(`TB_REFILE, bre);
        $readmemh(`TB_IMFILE, bim);
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
        $display("B_FFT %0d cycles and B_ENG %0d cycles for %0d frames (%0d and %0d per frame)",
                 c_fft, c_eng, `BSP_S, c_fft / `BSP_S, c_eng / `BSP_S);
        for (i = 0; i < 8; i = i + 1)  if (occ_top[i] != 0) $display("OCC,bsp_top,%0d,%0d", i, occ_top[i]);
        for (i = 0; i < 16; i = i + 1) if (occ_eng[i] != 0) $display("OCC,bsp_engine,%0d,%0d", i, occ_eng[i]);
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
