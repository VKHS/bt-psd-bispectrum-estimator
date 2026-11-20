// fd6_fir.v
// 6-tap fractional-delay FIR filter (real-valued).
// Coefficients correspond to h(n) in eq. (6):
//   {0.125, -0.2122, 0.6366, 0.6366, -0.2122, 0.125}
//
// This is a simple direct-form FIR (not folded). To match the paper's
// folded implementation, you'd share multipliers over time and/or
// reuse CRU. :contentReference[oaicite:7]{index=7}

`ifndef FD6_FIR_V
`define FD6_FIR_V

module fd6_fir #(
    parameter DATA_W = 16,
    parameter FRAC_W = 14
)(
    input  wire                     clk,
    input  wire                     rst_n,
    input  wire                     in_valid,
    input  wire signed [DATA_W-1:0] in_sample,
    output reg                      out_valid,
    output reg  signed [DATA_W-1:0] out_sample
);
    // Fixed-point coefficients (Q(FRAC_W))
    localparam signed [DATA_W-1:0] C0 = $rtoi( 0.1250  * (1<<FRAC_W));
    localparam signed [DATA_W-1:0] C1 = $rtoi(-0.2122 * (1<<FRAC_W));
    localparam signed [DATA_W-1:0] C2 = $rtoi( 0.6366 * (1<<FRAC_W));
    localparam signed [DATA_W-1:0] C3 = $rtoi( 0.6366 * (1<<FRAC_W));
    localparam signed [DATA_W-1:0] C4 = $rtoi(-0.2122 * (1<<FRAC_W));
    localparam signed [DATA_W-1:0] C5 = $rtoi( 0.1250 * (1<<FRAC_W));

    reg signed [DATA_W-1:0] x[0:5];

    integer i;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (i=0; i<6; i=i+1)
                x[i] <= '0;
            out_sample <= '0;
            out_valid  <= 1'b0;
        end else begin
            if (in_valid) begin
                // Shift register
                x[5] <= x[4];
                x[4] <= x[3];
                x[3] <= x[2];
                x[2] <= x[1];
                x[1] <= x[0];
                x[0] <= in_sample;

                // FIR dot product
                // Sum width: DATA_W + FRAC_W + log2(6) is safe
                // Here we compute combinationally and register at output.
                out_sample <=
                    (C0 * x[0] >>> FRAC_W) +
                    (C1 * x[1] >>> FRAC_W) +
                    (C2 * x[2] >>> FRAC_W) +
                    (C3 * x[3] >>> FRAC_W) +
                    (C4 * x[4] >>> FRAC_W) +
                    (C5 * x[5] >>> FRAC_W);

                out_valid <= 1'b1;
            end else begin
                out_valid <= 1'b0;
            end
        end
    end
endmodule

`endif
