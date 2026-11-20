// bt_estimator_top.sv
//
// Configurable Blackman–Tukey PSD / bispectrum estimator top.
// This is where we make "how necessary is this many?" concrete:
//
//   - If ENABLE_FD = 0  -> FD filter is bypassed
//   - If ENABLE_WINDOW = 0 -> window is bypassed (rectangular)
//   - If ENABLE_BISPECTRUM = 0 -> no bispectrum engine
//
// Build *one* design and compile it in different profiles. 

`timescale 1ns/1ps

import bt_config_pkg::*;

module bt_estimator_top #(
    parameter int W         = BT_WL,
    parameter int FRAC_BITS = BT_FRAC_BITS,
    parameter int LOGN      = BT_LOGN,
    parameter bit P_ENABLE_FD          = ENABLE_FD,
    parameter bit P_ENABLE_WINDOW      = ENABLE_WINDOW,
    parameter bit P_ENABLE_BISPECTRUM  = ENABLE_BISPECTRUM
)(
    input  wire                     clk,
    input  wire                     rst_n,

    // Streaming real input samples
    input  wire                     in_valid,
    input  wire signed [W-1:0]      in_sample,
    input  wire                     in_last,    // 1 on last sample of a frame of N=2^LOGN

    // Control to start PSD estimation on loaded frame
    input  wire                     start_frame,   // pulse when one frame is ready

    // PSD output: N bins of |X[k]|^2, sequential
    input  wire                     psd_read_en,
    output reg                      psd_valid,
    output reg  [W-1:0]             psd_bin,
    output reg                      psd_last,

    // Optional: simple bispectrum output interface (very simplified)
    input  wire                     bispec_req,    // request one bispectrum triple evaluation
    input  wire [LOGN-1:0]          f1_idx,
    input  wire [LOGN-1:0]          f2_idx,
    output reg                      bispec_valid,
    output reg  signed [W-1:0]      bispec_re,
    output reg  signed [W-1:0]      bispec_im
);

    localparam int N = (1 << LOGN);

    // ============================================================
    // 1) FRACTIONAL DELAY (optional)
    // ============================================================

    wire                 fd_valid;
    wire signed [W-1:0]  fd_sample;

    generate
        if (P_ENABLE_FD) begin : GEN_FD
            // Use the 6-tap bidirectional FD filter (eq. (6)). 
            fd6_fir #(
                .DATA_W (W),
                .FRAC_W (FRAC_BITS)
            ) u_fd6 (
                .clk        (clk),
                .rst_n      (rst_n),
                .in_valid   (in_valid),
                .in_sample  (in_sample),
                .out_valid  (fd_valid),
                .out_sample (fd_sample)
            );
        end else begin : GEN_FD_BYPASS
            assign fd_valid  = in_valid;
            assign fd_sample = in_sample;
        end
    endgenerate

    // ============================================================
    // 2) WINDOWING (optional)
    // ============================================================

    reg  [LOGN-1:0] win_index;
    wire            win_valid;
    wire signed [W-1:0] win_sample;

    // index inside frame for window ROM; increment when fd_valid
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            win_index <= '0;
        end else if (fd_valid) begin
            if (in_last) begin
                win_index <= '0;
            end else begin
                win_index <= win_index + 1'b1;
            end
        end
    end

    generate
        if (P_ENABLE_WINDOW) begin : GEN_WIN
            // e.g. time-domain Hann / Hamming / Blackman window
            // (coeffs in ROM generated from Table 4 selection). 
            window_fir #(
                .DATA_W   (W),
                .FRAC_W   (FRAC_BITS),
                .ADDR_W   (LOGN),
                .WIN_FILE ("window_hann.mem") // or hamming/blackman
            ) u_win (
                .clk        (clk),
                .rst_n      (rst_n),
                .in_valid   (fd_valid),
                .in_index   (win_index),
                .in_sample  (fd_sample),
                .out_valid  (win_valid),
                .out_sample (win_sample)
            );
        end else begin : GEN_WIN_BYPASS
            assign win_valid  = fd_valid;
            assign win_sample = fd_sample;
        end
    endgenerate

    // ============================================================
    // 3) FFT CORE (memory-based, safe scaling)
    // ============================================================

    reg                  fft_load_valid;
    reg  signed [W-1:0]  fft_load_re, fft_load_im;
    reg                  fft_load_last;
    reg                  fft_start;

    wire                 fft_busy;
    wire                 fft_done;
    wire                 fft_out_valid;
    wire signed [W-1:0]  fft_out_re, fft_out_im;

    fft_core_iterative #(
        .W         (W),
        .FRAC_BITS (FRAC_BITS),
        .LOGN      (LOGN)
    ) u_fft (
        .clk        (clk),
        .rst_n      (rst_n),
        .load_valid (fft_load_valid),
        .load_re    (fft_load_re),
        .load_im    (fft_load_im),
        .load_last  (fft_load_last),
        .start_fft  (fft_start),
        .busy       (fft_busy),
        .done       (fft_done),
        .read_en    (psd_read_en),
        .read_valid (fft_out_valid),
        .out_re     (fft_out_re),
        .out_im     (fft_out_im)
    );

    // Simple frame loader FSM from windowed samples into FFT
    localparam P_IDLE  = 2'd0;
    localparam P_LOAD  = 2'd1;
    localparam P_WAIT  = 2'd2;

    reg [1:0]           p_state;
    reg [LOGN-1:0]      load_count;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            p_state        <= P_IDLE;
            fft_load_valid <= 1'b0;
            fft_load_re    <= '0;
            fft_load_im    <= '0;
            fft_load_last  <= 1'b0;
            fft_start      <= 1'b0;
            load_count     <= '0;
        end else begin
            fft_load_valid <= 1'b0;
            fft_load_last  <= 1'b0;
            fft_start      <= 1'b0;

            case (p_state)
                P_IDLE: begin
                    load_count <= 0;
                    if (start_frame) begin
                        p_state <= P_LOAD;
                    end
                end

                P_LOAD: begin
                    if (win_valid) begin
                        fft_load_valid <= 1'b1;
                        fft_load_re    <= win_sample;
                        fft_load_im    <= '0;
                        fft_load_last  <= (load_count == N-1);
                        load_count     <= load_count + 1'b1;

                        if (load_count == N-1) begin
                            p_state   <= P_WAIT;
                            fft_start <= 1'b1;
                        end
                    end
                end

                P_WAIT: begin
                    if (fft_done) begin
                        p_state <= P_IDLE;
                    end
                end

                default: p_state <= P_IDLE;
            endcase
        end
    end

    // ============================================================
    // 4) |X|^2 for PSD
    // ============================================================

    // Here we use a simple multiplier-based squarer; if you want to
    // match the sum-of-odds method of eq. (26), replace this with
    // mag_sq_odd. 
    wire [2*W-1:0] re2 = fft_out_re * fft_out_re;
    wire [2*W-1:0] im2 = fft_out_im * fft_out_im;
    wire [2*W  :0] mag2_ext = re2 + im2;

    // Truncate/scale down to W bits (you can refine scaling)
    wire [W-1:0] mag2_q = mag2_ext[W-1:0];

    reg [LOGN-1:0] psd_idx;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            psd_valid <= 1'b0;
            psd_bin   <= '0;
            psd_last  <= 1'b0;
            psd_idx   <= '0;
        end else begin
            psd_valid <= 1'b0;
            psd_last  <= 1'b0;

            if (fft_out_valid && psd_read_en) begin
                psd_bin   <= mag2_q;
                psd_valid <= 1'b1;

                if (psd_idx == N-1) begin
                    psd_last <= 1'b1;
                    psd_idx  <= 0;
                end else begin
                    psd_idx  <= psd_idx + 1'b1;
                end
            end
        end
    end

    // ============================================================
    // 5) (Optional) Bispectrum engine
    // ============================================================

    // For simplicity, this block:
    //   - uses the *current* FFT frame in RAM (fft_core_iterative)
    //   - when bispec_req is asserted, reads F(f1), F(f2), F(f1+f2)
    //   - runs one triple product through bispectrum_engine
    //
    // In real hardware you'd average bispectrum values over many
    // segments, similar to the PSD Save/Add unit discussed in eq. (27). 

    // We assume fft_core_iterative internally stores its frame in RAM
    // and lets us re-read bins by address. If not, you'd add a RAM
    // wrapper for that.

    generate
        if (P_ENABLE_BISPECTRUM) begin : GEN_BISPEC

            // A *very* simplified view: we treat fft_core as if we can
            // access stored bins by index through some side-port.
            // For now we'll model that as wires; in your concrete design
            // you would map these to a RAM or to a retained buffer.

            // Placeholder "readback" interface – you should connect
            // this to your FFT data RAM.
            function automatic signed [W-1:0] fft_re_at;
                input [LOGN-1:0] idx;
                begin
                    // TODO: connect to actual RAM / storage.
                    fft_re_at = '0;
                end
            endfunction

            function automatic signed [W-1:0] fft_im_at;
                input [LOGN-1:0] idx;
                begin
                    fft_im_at = '0;
                end
            endfunction

            wire signed [W-1:0] F1_re = fft_re_at(f1_idx);
            wire signed [W-1:0] F1_im = fft_im_at(f1_idx);
            wire signed [W-1:0] F2_re = fft_re_at(f2_idx);
            wire signed [W-1:0] F2_im = fft_im_at(f2_idx);
            wire signed [W-1:0] F3_re = fft_re_at((f1_idx + f2_idx) & (N-1));
            wire signed [W-1:0] F3_im = fft_im_at((f1_idx + f2_idx) & (N-1));

            wire signed [W-1:0] B_re_w, B_im_w;
            wire                 B_valid_w;

            bispectrum_engine #(
                .DATA_W (W),
                .FRAC_W (FRAC_BITS)
            ) u_bis (
                .clk       (clk),
                .rst_n     (rst_n),
                .in_valid  (bispec_req),
                .f1_r      (F1_re),
                .f1_i      (F1_im),
                .f2_r      (F2_re),
                .f2_i      (F2_im),
                .f3_r      (F3_re),
                .f3_i      (F3_im),
                .out_valid (B_valid_w),
                .b_r       (B_re_w),
                .b_i       (B_im_w)
            );

            always @(posedge clk or negedge rst_n) begin
                if (!rst_n) begin
                    bispec_valid <= 1'b0;
                    bispec_re    <= '0;
                    bispec_im    <= '0;
                end else begin
                    bispec_valid <= 1'b0;
                    if (B_valid_w) begin
                        bispec_valid <= 1'b1;
                        bispec_re    <= B_re_w;
                        bispec_im    <= B_im_w;
                    end
                end
            end

        end else begin : GEN_BISPEC_NONE
            always @(posedge clk or negedge rst_n) begin
                if (!rst_n) begin
                    bispec_valid <= 1'b0;
                    bispec_re    <= '0;
                    bispec_im    <= '0;
                end else begin
                    bispec_valid <= 1'b0;
                    bispec_re    <= '0;
                    bispec_im    <= '0;
                end
            end
        end
    endgenerate

endmodule
