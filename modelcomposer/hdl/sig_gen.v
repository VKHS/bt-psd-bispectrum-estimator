// Test-signal generator: a sine tone (phase accumulator and a 1024-entry sine table) plus pseudo-random noise (32-bit LFSR),
// summed with saturation.  One new sample per clock cycle in which en = 1; every register holds while en = 0.
//
//   sample n:   x = sat( ROM[ phase_n >> (PHW-10) ]  +  (signed top WL bits of lfsr_n) >>> noise_shift )
//   phase_n  = n * pinc  (mod 2^PHW),   lfsr_n = state after n updates from SEED
//
// It plays the role that the Simulink Sine Wave and Random Number sources, or the AMD DDS Compiler, play in a Model Composer
// design, as plain Verilog: it is NOT the AMD DDS Compiler IP.  The output of sample 0 is visible two enabled clock edges after
// reset (valid = 1), then a new sample at every enabled edge.  Bit-exact with btpsd.siggen.generate (tb/tb_sig_gen.v).
//
//   pinc          phase increment, tone frequency = pinc / 2^PHW cycles per sample (btpsd.siggen.pinc_for(bin) for bin / 1024)
//   noise_en      0: tone only
//   noise_shift   noise amplitude 2^(WL-1-noise_shift), uniform, 0 .. WL-1
//   rst_n         synchronous, active low; the registers also start in the reset state without it
`timescale 1ns / 1ps

module sig_gen #(
    parameter WL   = 12,
    parameter PHW  = 32,
    parameter [31:0] SEED = 32'hA5A51F3D
) (
    input                      clk,
    input                      rst_n,
    input                      en,
    input      [PHW-1:0]       pinc,
    input                      noise_en,
    input      [3:0]           noise_shift,
    output signed [WL-1:0]     x,
    output                     valid
);
    reg [PHW-1:0] phase = 0;
    reg [31:0]    lfsr  = SEED;
    reg           v1    = 1'b0;
    reg signed [WL-1:0] nz_q = 0;
    reg signed [WL-1:0] x_r  = 0;
    reg                 v2   = 1'b0;
    assign x     = x_r;
    assign valid = v2;

    wire [WL-1:0] rom_q;
    sig_rom #(.WL(WL)) u_rom (.clk(clk), .en(en), .addr(phase[PHW-1:PHW-10]), .dout(rom_q));

    wire signed [WL-1:0] nz_raw = lfsr[31:32-WL];
    wire signed [WL-1:0] nz_s   = nz_raw >>> noise_shift;
    wire signed [WL-1:0] tone   = rom_q;
    wire signed [WL:0]   sum    = tone + nz_q;
    localparam signed [WL:0] MAXV = (1 << (WL - 1)) - 1;
    localparam signed [WL:0] MINV = -(1 << (WL - 1));
    localparam        [WL-1:0] MAXW = (1 << (WL - 1)) - 1;      // the same limits as WL-bit patterns
    localparam        [WL-1:0] MINW = (1 << (WL - 1));

    always @(posedge clk) begin
        if (!rst_n) begin
            phase <= 0;  lfsr <= SEED;  v1 <= 1'b0;  nz_q <= 0;  x_r <= 0;  v2 <= 1'b0;
        end else if (en) begin
            // stage 1: table read (inside sig_rom), noise word, advance phase and LFSR
            nz_q  <= noise_en ? nz_s : {WL{1'b0}};
            phase <= phase + pinc;
            lfsr  <= (lfsr >> 1) ^ (lfsr[0] ? 32'h80200003 : 32'h0);
            v1    <= 1'b1;
            // stage 2: sum and saturate
            x_r   <= (sum > MAXV) ? MAXW : (sum < MINV) ? MINW : sum[WL-1:0];
            v2    <= v1;
        end
    end
endmodule
