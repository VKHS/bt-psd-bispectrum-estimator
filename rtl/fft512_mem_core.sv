// fft512_mem_core.sv  (revised)
//
// Memory-based 512-point iterative radix-2 FFT core with:
//   - dual-port RAM for real/imag
//   - proper Cooley–Tukey DIT schedule
//   - per-stage safe scaling (>>1) so total gain = 1/N (eq. (12)) :contentReference[oaicite:2]{index=2}
//   - twiddle_rom providing W_N^k
//
// One butterfly per clock (not counting load/output phases).
// You can wrap this in your BT estimator top-level from before.
//
// NOTE:
//   - For clarity this uses multipliers inside the butterfly instead of a
//     CORDIC-II CRU; structurally it matches the paper's memory-based
//     FFT with a central rotator, but not multiplierless yet. :contentReference[oaicite:3]{index=3}

`include "dual_port_ram.sv"
`include "radix2_butterfly.sv"
`include "twiddle_rom.sv"

module fft512_mem_core #(
    parameter int N_POINTS  = 512,
    parameter int LOGN      = 9,    // log2(512)
    parameter int DATA_W    = 16,
    parameter int ADDR_W    = LOGN,
    parameter int TWID_W    = 16
) (
    input  logic                      clk,
    input  logic                      rst_n,

    // Streamed real input samples, imag = 0
    input  logic                      in_valid,
    output logic                      in_ready,
    input  logic signed [DATA_W-1:0]  in_sample,

    // Streamed complex output samples
    output logic                      out_valid,
    input  logic                      out_ready,
    output logic signed [DATA_W-1:0]  out_re,
    output logic signed [DATA_W-1:0]  out_im,

    // Control
    input  logic                      start,
    output logic                      done
);

    // ----------------------------------------------------------------
    // RAMs for real and imag parts
    // ----------------------------------------------------------------
    logic                     we_a, we_b;
    logic [ADDR_W-1:0]        addr_a, addr_b;
    logic signed [DATA_W-1:0] din_re_a, din_re_b;
    logic signed [DATA_W-1:0] dout_re_a, dout_re_b;
    logic signed [DATA_W-1:0] din_im_a, din_im_b;
    logic signed [DATA_W-1:0] dout_im_a, dout_im_b;

    dual_port_ram #(
        .ADDR_W(ADDR_W),
        .DATA_W(DATA_W)
    ) ram_real (
        .clk   (clk),
        .we_a  (we_a),
        .addr_a(addr_a),
        .din_a (din_re_a),
        .dout_a(dout_re_a),
        .we_b  (we_b),
        .addr_b(addr_b),
        .din_b (din_re_b),
        .dout_b(dout_re_b)
    );

    dual_port_ram #(
        .ADDR_W(ADDR_W),
        .DATA_W(DATA_W)
    ) ram_imag (
        .clk   (clk),
        .we_a  (we_a),
        .addr_a(addr_a),
        .din_a (din_im_a),
        .dout_a(dout_im_a),
        .we_b  (we_b),
        .addr_b(addr_b),
        .din_b (din_im_b),
        .dout_b(dout_im_b)
    );

    // ----------------------------------------------------------------
    // FSM states
    // ----------------------------------------------------------------
    typedef enum logic [1:0] {
        S_IDLE,
        S_LOAD,
        S_FFT,
        S_OUTPUT
    } state_t;

    state_t state, state_n;

    // Counters
    logic [ADDR_W-1:0] load_idx, load_idx_n;
    logic [ADDR_W-1:0] out_idx,  out_idx_n;

    // FFT stage / group / j counters
    logic [LOGN-1:0]   stage, stage_n;
    logic [ADDR_W-1:0] group, group_n;  // 0 .. N/m - 1
    logic [ADDR_W-1:0] j,     j_n;      // 0 .. half_m - 1

    // ----------------------------------------------------------------
    // Twiddle ROM
    // ----------------------------------------------------------------
    localparam int TWID_ADDR_W = LOGN - 1; // log2(N/2) = 8
    logic [TWID_ADDR_W-1:0]    tw_addr;
    logic signed [TWID_W-1:0]  tw_re_full, tw_im_full;
    logic signed [DATA_W-1:0]  tw_re, tw_im;

    twiddle_rom #(
        .N      (N_POINTS),
        .DATA_W (TWID_W),
        .ADDR_W (TWID_ADDR_W),
        .RE_FILE("twiddle_re_512_q15.mem"),
        .IM_FILE("twiddle_im_512_q15.mem")
    ) u_tw (
        .addr  (tw_addr),
        .tw_re (tw_re_full),
        .tw_im (tw_im_full)
    );

    // Truncate / align twiddle to DATA_W
    assign tw_re = tw_re_full[TWID_W-1 -: DATA_W];
    assign tw_im = tw_im_full[TWID_W-1 -: DATA_W];

    // ----------------------------------------------------------------
    // Butterfly
    // ----------------------------------------------------------------
    logic signed [DATA_W-1:0] a_re, a_im, b_re, b_im;
    logic signed [DATA_W-1:0] y0_re, y0_im, y1_re, y1_im;

    radix2_butterfly #(
        .DATA_W(DATA_W)
    ) u_bfly (
        .a_re (a_re),
        .a_im (a_im),
        .b_re (b_re),
        .b_im (b_im),
        .w_re (tw_re),
        .w_im (tw_im),
        .y0_re(y0_re),
        .y0_im(y0_im),
        .y1_re(y1_re),
        .y1_im(y1_im)
    );

    // Safe scaling: divide outputs by 2 every stage (right shift 1).
    logic signed [DATA_W-1:0] y0_re_sc, y0_im_sc, y1_re_sc, y1_im_sc;
    assign y0_re_sc = y0_re >>> 1;
    assign y0_im_sc = y0_im >>> 1;
    assign y1_re_sc = y1_re >>> 1;
    assign y1_im_sc = y1_im >>> 1;

    // ----------------------------------------------------------------
    // State and counter registers
    // ----------------------------------------------------------------
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state     <= S_IDLE;
            load_idx  <= '0;
            out_idx   <= '0;
            stage     <= '0;
            group     <= '0;
            j         <= '0;
        end else begin
            state     <= state_n;
            load_idx  <= load_idx_n;
            out_idx   <= out_idx_n;
            stage     <= stage_n;
            group     <= group_n;
            j         <= j_n;
        end
    end

    // ----------------------------------------------------------------
    // Combinational control
    // ----------------------------------------------------------------
    always_comb begin
        // defaults
        state_n     = state;
        load_idx_n  = load_idx;
        out_idx_n   = out_idx;
        stage_n     = stage;
        group_n     = group;
        j_n         = j;

        in_ready    = 1'b0;
        out_valid   = 1'b0;
        out_re      = '0;
        out_im      = '0;
        done        = 1'b0;

        we_a        = 1'b0;
        we_b        = 1'b0;
        addr_a      = '0;
        addr_b      = '0;
        din_re_a    = '0;
        din_im_a    = '0;
        din_re_b    = '0;
        din_im_b    = '0;

        // Derived for FFT
        logic [LOGN:0]  m;
        logic [LOGN:0]  half_m;
        logic [LOGN:0]  groups;
        logic [ADDR_W-1:0] addr_a_c, addr_b_c;
        logic [TWID_ADDR_W-1:0] tw_addr_c;
        logic [ADDR_W-1:0] k_step;

        m       = (1 << (stage + 1));      // 2^(stage+1)
        half_m  = (1 << stage);            // 2^stage
        groups  = (N_POINTS >> (stage+1)); // N / m
        k_step  = (N_POINTS >> (stage+1)); // same as groups

        // Twiddle index (0..N/2-1)
        tw_addr_c = j * k_step;
        tw_addr   = tw_addr_c;

        // Connect RAM outputs to butterfly inputs
        a_re = dout_re_a;
        a_im = dout_im_a;
        b_re = dout_re_b;
        b_im = dout_im_b;

        case (state)
            // ----------------------------------------
            S_IDLE: begin
                if (start) begin
                    state_n    = S_LOAD;
                    load_idx_n = '0;
                end
            end

            // ----------------------------------------
            S_LOAD: begin
                in_ready = 1'b1;
                if (in_valid) begin
                    addr_a    = load_idx;
                    din_re_a  = in_sample;
                    din_im_a  = '0;
                    we_a      = 1'b1;

                    if (load_idx == N_POINTS-1) begin
                        // transition to FFT
                        state_n    = S_FFT;
                        stage_n    = '0;
                        group_n    = '0;
                        j_n        = '0;
                    end else begin
                        load_idx_n = load_idx + 1;
                    end
                end
            end

            // ----------------------------------------
            S_FFT: begin
                // Addresses for this butterfly: a and b
                addr_a_c = group * m + j;
                addr_b_c = addr_a_c + half_m[ADDR_W-1:0];

                addr_a   = addr_a_c;
                addr_b   = addr_b_c;

                // Write back scaled butterfly outputs
                we_a     = 1'b1;
                we_b     = 1'b1;
                din_re_a = y0_re_sc;
                din_im_a = y0_im_sc;
                din_re_b = y1_re_sc;
                din_im_b = y1_im_sc;

                // End-of-loop flags
                logic last_j, last_group, last_stage, last_bfly_stage;
                last_j          = (j == half_m[ADDR_W-1:0] - 1);
                last_group      = (group == groups[ADDR_W-1:0] - 1);
                last_stage      = (stage == LOGN-1);
                last_bfly_stage = last_j && last_group;

                // Counter updates
                if (last_bfly_stage) begin
                    if (last_stage) begin
                        // FFT all stages done
                        state_n   = S_OUTPUT;
                        out_idx_n = '0;
                    end else begin
                        // Move to next stage
                        stage_n = stage + 1;
                        group_n = '0;
                        j_n     = '0;
                    end
                end else begin
                    // Continue within this stage
                    if (last_j) begin
                        // next group
                        group_n = group + 1;
                        j_n     = '0;
                    end else begin
                        // next j in same group
                        j_n = j + 1;
                    end
                end
            end

            // ----------------------------------------
            S_OUTPUT: begin
                if (out_ready) begin
                    addr_a   = out_idx;
                    out_re   = dout_re_a;
                    out_im   = dout_im_a;
                    out_valid= 1'b1;

                    if (out_idx == N_POINTS-1) begin
                        state_n = S_IDLE;
                        done    = 1'b1;
                    end else begin
                        out_idx_n = out_idx + 1;
                    end
                end
            end

            default: state_n = S_IDLE;
        endcase
    end

endmodule
