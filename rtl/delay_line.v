// N-stage delay line of a W-bit word (N >= 1).
module delay_line #(
    parameter W = 8,
    parameter N = 1
) (
    input          clk,
    input  [W-1:0] d,
    output [W-1:0] q
);
    reg [W-1:0] sr [0:N-1];
    integer i;
    always @(posedge clk) begin
        sr[0] <= d;
        for (i = 1; i < N; i = i + 1) sr[i] <= sr[i-1];
    end
    assign q = sr[N-1];
endmodule
