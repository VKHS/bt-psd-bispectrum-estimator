// Simple dual-port RAM: one synchronous read port and one write port.
// A read of the address written in the same cycle returns the old word (read first).
module ram_1r1w #(
    parameter DW = 32,
    parameter AW = 6
) (
    input               clk,
    input               we,
    input  [AW-1:0]     waddr,
    input  [DW-1:0]     wdata,
    input  [AW-1:0]     raddr,
    output reg [DW-1:0] rdata
);
    reg [DW-1:0] mem [0:(1<<AW)-1];
    always @(posedge clk) begin
        if (we) mem[waddr] <= wdata;
        rdata <= mem[raddr];
    end
endmodule
