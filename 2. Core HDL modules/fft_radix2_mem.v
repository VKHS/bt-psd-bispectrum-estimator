// fft_radix2_mem.v
// Sequential in-place radix-2 FFT with internal memory,
// using cru_shared for twiddle multiplication.
//
// This is a "core" FFT engine; the Blackman–Tukey estimator
// will sit on top of it and feed/consume samples in blocks.

`ifndef FFT_RADIX2_MEM_V
`define FFT_RADIX2_MEM_V

`include "complex_fixed.v"

module fft_radix2_mem #(
    parameter DATA_W  = 16,
    parameter FRAC_W  = 14,
    parameter N_LOG2  = 9   // log2(N); N=2^N_LOG2 (use 9 for 512)
)(
    input  wire                     clk,
    input  wire                     rst_n,

    // Input block interface: N samples, one per cycle when in_valid && in_ready
    input  wire                     in_valid,
    output wire                     in_ready,
    input  wire signed [DATA_W-1:0] in_r,
    input  wire signed [DATA_W-1:0] in_i,

    // Output block interface: N FFT results, one per cycle
    output reg                      out_valid,
    input  wire                     out_ready,
    output reg  signed [DATA_W-1:0] out_r,
    output reg  signed [DATA_W-1:0] out_i,

    output reg                      done   // pulses high for 1 cycle after last output
);

    localparam N      = (1 << N_LOG2);
    localparam ST_IDLE   = 3'd0;
    localparam ST_LOAD   = 3'd1;
    localparam ST_STAGE  = 3'd2;
    localparam ST_OUTPUT = 3'd3;

    reg [2:0] state, next_state;

    // Internal memory for complex samples
    reg signed [DATA_W-1:0] mem_r [0:N-1];
    reg signed [DATA_W-1:0] mem_i [0:N-1];

    // Counters for loading, stages, butterflies, and output
    reg [N_LOG2-1:0] load_idx;
    reg [N_LOG2-1:0] out_idx;

    reg [N_LOG2-1:0] stage;         // 0..N_LOG2-1
    reg [N_LOG2-1:0] butterfly;     // 0..N/2-1
    reg [N_LOG2-1:0] m;             // size of sub-FFT at this stage
    reg [N_LOG2-1:0] half_m;
    reg [N_LOG2-1:0] j;
    reg [N_LOG2-1:0] k;

    // CRU connection
    wire         cru_req_ready;
    reg          cru_req_valid;
    reg signed [DATA_W-1:0] cru_xr, cru_xi;
    reg [N_LOG2-1:0]        cru_angle_idx;
    wire         cru_resp_valid;
    reg          cru_resp_ready;
    wire signed [DATA_W-1:0] cru_yr, cru_yi;

    // Instantiate CRU
    cru_shared #(
        .DATA_W  (DATA_W),
        .FRAC_W  (FRAC_W),
        .ANGLE_W (N_LOG2),
        .N_TW    (N)
    ) u_cru (
        .clk        (clk),
        .rst_n      (rst_n),
        .req_valid  (cru_req_valid),
        .req_ready  (cru_req_ready),
        .xr         (cru_xr),
        .xi         (cru_xi),
        .angle_idx  (cru_angle_idx),
        .resp_valid (cru_resp_valid),
        .resp_ready (cru_resp_ready),
        .yr         (cru_yr),
        .yi         (cru_yi)
    );

    // Input ready when loading and not finished
    assign in_ready = (state == ST_LOAD);

    // General state machine
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state     <= ST_IDLE;
            load_idx  <= 0;
            stage     <= 0;
            butterfly <= 0;
            m         <= 0;
            half_m    <= 0;
            j         <= 0;
            k         <= 0;
            out_idx   <= 0;
            out_valid <= 1'b0;
            done      <= 1'b0;
            cru_req_valid  <= 1'b0;
            cru_resp_ready <= 1'b0;
        end else begin
            state <= next_state;
            done  <= 1'b0;

            case (state)
                ST_IDLE: begin
                    // Wait for first input sample to start loading
                    if (in_valid) begin
                        mem_r[0] <= in_r;
                        mem_i[0] <= in_i;
                        load_idx <= 1;
                    end
                end

                ST_LOAD: begin
                    if (in_valid && in_ready) begin
                        mem_r[load_idx] <= in_r;
                        mem_i[load_idx] <= in_i;
                        load_idx <= load_idx + 1'b1;
                    end
                end

                ST_STAGE: begin
                    // Sequentially handle butterflies
                    // When we start a stage, m and half_m are set.

                    // If no active CRU request, launch one
                    if (!cru_req_valid && !cru_resp_valid) begin
                        // Butterfly indices: positions k+j and k+j+half_m
                        // Twiddle angle index: (j * N / m)
                        // Access mem combinationally and send v to CRU
                        cru_xr        <= mem_r[k + j + half_m];
                        cru_xi        <= mem_i[k + j + half_m];
                        cru_angle_idx <= (j << (N_LOG2 - $clog2(m)));
                        cru_req_valid <= 1'b1;
                        cru_resp_ready <= 1'b1;
                    end else if (cru_req_valid && cru_req_ready) begin
                        // Request accepted
                        cru_req_valid <= 1'b0;
                    end

                    // When CRU responds, compute butterfly and write back
                    if (cru_resp_valid && cru_resp_ready) begin
                        // u = mem[k+j], t = cru_yr/cru_yi
                        // y0 = u + t, y1 = u - t
                        // We ignore safe scaling here – can be added as >>1 each stage.
                        mem_r[k+j]          <= mem_r[k+j] + cru_yr;
                        mem_i[k+j]          <= mem_i[k+j] + cru_yi;
                        mem_r[k+j+half_m]   <= mem_r[k+j] - cru_yr;
                        mem_i[k+j+half_m]   <= mem_i[k+j] - cru_yi;

                        // Advance butterfly indices
                        if (j == half_m - 1) begin
                            j <= 0;
                            if (k == N - m) begin
                                // Stage finished
                                butterfly <= 0;
                            end else begin
                                k <= k + m;
                            end
                        end else begin
                            j <= j + 1'b1;
                        end
                    end
                end

                ST_OUTPUT: begin
                    if (!out_valid && out_idx < N) begin
                        out_r    <= mem_r[out_idx];
                        out_i    <= mem_i[out_idx];
                        out_valid <= 1'b1;
                    end else if (out_valid && out_ready) begin
                        out_idx  <= out_idx + 1'b1;
                        out_valid <= 1'b0;
                        if (out_idx == N-1) begin
                            done <= 1'b1;
                        end
                    end
                end

                default: ;
            endcase
        end
    end

    // Next state logic and stage management
    always @(*) begin
        next_state = state;

        case (state)
            ST_IDLE: begin
                if (in_valid)
                    next_state = ST_LOAD;
            end

            ST_LOAD: begin
                if (load_idx == N-1 && in_valid && in_ready) begin
                    // Loaded all N samples
                    next_state = ST_STAGE;
                end
            end

            ST_STAGE: begin
                // When we've finished all stages, move to output
                if (stage == N_LOG2 && butterfly == 0 && j == 0) begin
                    next_state = ST_OUTPUT;
                end
            end

            ST_OUTPUT: begin
                if (done)
                    next_state = ST_IDLE;
            end

            default: next_state = ST_IDLE;
        endcase
    end

    // Stage control (simple version – you can refine for correctness)
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            stage     <= 0;
            m         <= 2;
            half_m    <= 1;
            k         <= 0;
            j         <= 0;
        end else if (state == ST_STAGE) begin
            // This is skeletal; in a full implementation you’d
            // advance stage when butterflies for current stage are done.
            if (cru_resp_valid && cru_resp_ready && k == (N - m) && j == (half_m - 1)) begin
                // Move to next stage
                stage  <= stage + 1'b1;
                m      <= m << 1;
                half_m <= half_m << 1;
                k      <= 0;
                j      <= 0;
            end
        end else if (state == ST_LOAD) begin
            // Re-init stage values when we start FFT
            stage  <= 1;
            m      <= 2;
            half_m <= 1;
            k      <= 0;
            j      <= 0;
        end
    end

endmodule

`endif
