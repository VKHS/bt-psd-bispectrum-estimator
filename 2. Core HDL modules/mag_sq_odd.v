// mag_sq_odd.v
// Compute |X|^2 using (Re^2 + Im^2), where each square is realized
// by sum of odd numbers: n^2 = sum_{i=1..n}(2i-1). :contentReference[oaicite:10]{index=10}
//
// For clarity we treat the magnitude inputs as non-negative integers
// (take absolute value before feeding).

`ifndef MAG_SQ_ODD_V
`define MAG_SQ_ODD_V

module mag_sq_odd #(
    parameter DATA_W = 16
)(
    input  wire                 clk,
    input  wire                 rst_n,
    input  wire                 in_valid,
    input  wire signed [DATA_W-1:0] in_r,
    input  wire signed [DATA_W-1:0] in_i,
    output reg                  out_valid,
    output reg  [2*DATA_W-1:0]  out_mag_sq
);
    localparam ST_IDLE  = 2'd0;
    localparam ST_RSQ   = 2'd1;
    localparam ST_ISQ   = 2'd2;
    localparam ST_DONE  = 2'd3;

    reg [1:0] state;

    reg [DATA_W-1:0] n_r;      // |Re|
    reg [DATA_W-1:0] n_i;      // |Im|
    reg [DATA_W-1:0] cnt;
    reg [2*DATA_W-1:0] acc;
    reg [2*DATA_W-1:0] rsq, isq;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state      <= ST_IDLE;
            out_valid  <= 1'b0;
            out_mag_sq <= '0;
            n_r        <= '0;
            n_i        <= '0;
            cnt        <= '0;
            acc        <= '0;
            rsq        <= '0;
            isq        <= '0;
        end else begin
            out_valid <= 1'b0;

            case (state)
                ST_IDLE: begin
                    if (in_valid) begin
                        n_r <= in_r[DATA_W-1] ? -in_r : in_r;
                        n_i <= in_i[DATA_W-1] ? -in_i : in_i;
                        cnt <= 1;
                        acc <= 1; // first odd is 1
                        state <= ST_RSQ;
                    end
                end

                ST_RSQ: begin
                    if (cnt < n_r) begin
                        cnt <= cnt + 1'b1;
                        acc <= acc + (cnt<<1) - 1; // add next odd
                    end else begin
                        rsq <= acc;
                        // reset for imaginary
                        cnt <= 1;
                        acc <= 1;
                        state <= ST_ISQ;
                    end
                end

                ST_ISQ: begin
                    if (cnt < n_i) begin
                        cnt <= cnt + 1'b1;
                        acc <= acc + (cnt<<1) - 1;
                    end else begin
                        isq        <= acc;
                        out_mag_sq <= rsq + acc;
                        out_valid  <= 1'b1;
                        state      <= ST_DONE;
                    end
                end

                ST_DONE: begin
                    // One-cycle done state, then back to idle
                    state <= ST_IDLE;
                end

                default: state <= ST_IDLE;
            endcase
        end
    end
endmodule

`endif
