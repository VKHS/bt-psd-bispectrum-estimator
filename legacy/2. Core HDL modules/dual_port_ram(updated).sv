// dual_port_ram.sv  (updated)
//
// Dual-port RAM with synchronous write, asynchronous read.
// Simpler to use in the iterative FFT core.
//

module dual_port_ram #(
    parameter int ADDR_W = 9,
    parameter int DATA_W = 16
) (
    input  logic                  clk,
    // Port A
    input  logic                  we_a,
    input  logic [ADDR_W-1:0]     addr_a,
    input  logic [DATA_W-1:0]     din_a,
    output logic [DATA_W-1:0]     dout_a,
    // Port B
    input  logic                  we_b,
    input  logic [ADDR_W-1:0]     addr_b,
    input  logic [DATA_W-1:0]     din_b,
    output logic [DATA_W-1:0]     dout_b
);

    logic [DATA_W-1:0] mem [0:(1<<ADDR_W)-1];

    always_ff @(posedge clk) begin
        if (we_a) mem[addr_a] <= din_a;
        if (we_b) mem[addr_b] <= din_b;
    end

    // async read
    assign dout_a = mem[addr_a];
    assign dout_b = mem[addr_b];

endmodule
