// psd_accumulator.v
// Accumulate |X|^2 over multiple segments with right-shift scaling.
// You can use this for Blackman–Tukey averaging without explicit
// division: safe scaling inside FFT gives 1/L; this accumulator gives
// another 1/2^SCALE_SH right shift to approximate 1/K. :contentReference[oaicite:12]{index=12}

`ifndef PSD_ACCUMULATOR_V
`define PSD_ACCUMULATOR_V

module psd_accumulator #(
    parameter MAG_W      = 32,  // width of |X|^2
    parameter ACC_W      = 40,  // accumulator width
    parameter SCALE_SH   = 4    // shift right by SCALE_SH after accumulation
)(
    input  wire                 clk,
    input  wire                 rst_n,
    input  wire                 clear,      // clear accumulators
    input  wire                 add_valid,  // |X|^2 value valid
    input  wire [MAG_W-1:0]     mag_sq,
    output reg  [ACC_W-1:0]     acc_out,
    output reg                  acc_valid
);
    reg [ACC_W-1:0] acc_reg;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            acc_reg  <= '0;
            acc_out  <= '0;
            acc_valid <= 1'b0;
        end else begin
            if (clear) begin
                acc_reg <= '0;
                acc_valid <= 1'b0;
            end else if (add_valid) begin
                acc_reg <= acc_reg + mag_sq;
            end

            // For simplicity, we emit scaled value when clear is asserted
            if (clear) begin
                acc_out  <= acc_reg >> SCALE_SH;
                acc_valid <= 1'b1;
            end else begin
                acc_valid <= 1'b0;
            end
        end
    end
endmodule

`endif
