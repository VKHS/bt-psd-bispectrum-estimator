// Self-checking testbench of cru_rotate against btpsd.cordic.rotate_fixed (WL = 16, 9216 vectors).
// Run from the repository root:  make -C sim cru
`timescale 1ns / 1ps

module tb_cru_rotate;
    localparam WL = 16, N = 9216, LAT = 4;

    reg clk = 1'b0;
    always #5 clk = ~clk;

    reg [2*WL+9:0] vin  [0:N-1];
    reg [2*WL-1:0] vout [0:N-1];

    reg signed [WL-1:0] xr, xi;
    reg        [9:0]    q;
    wire signed [WL-1:0] yr, yi;

    cru_rotate #(.WL(WL), .ROUND(1)) dut (.clk(clk), .xr(xr), .xi(xi), .q(q), .yr(yr), .yi(yi));

    integer i, errs;
    initial begin
        $readmemh("sim/vectors/cru_wl16_in.mem", vin);
        $readmemh("sim/vectors/cru_wl16_out.mem", vout);
        errs = 0;
        xr = 0;  xi = 0;  q = 0;
        for (i = 0; i < N + LAT; i = i + 1) begin
            @(negedge clk);
            if (i >= LAT) begin
                if ({yr, yi} !== vout[i-LAT]) begin
                    errs = errs + 1;
                    if (errs <= 10) $display("MISMATCH vector %0d: got %h expected %h", i-LAT, {yr, yi}, vout[i-LAT]);
                end
            end
            if (i < N) {xr, xi, q} = vin[i];
        end
        if (errs == 0) $display("TB_RESULT: PASS tb_cru_rotate (%0d vectors)", N);
        else           $display("TB_RESULT: FAIL tb_cru_rotate (%0d mismatches)", errs);
        $finish;
    end
endmodule
