// bt_config_pkg.sv
//
// Configuration for the Blackman–Tukey PSD / bispectrum estimator.
//
// This answers "How necessary is this many?" by letting you *turn off*
// optional blocks (FD, window, bispectrum) and still build a working
// design. 
//
// PROFILES:
//   - Minimal PSD (no FD, no window, no bispectrum)
//   - PSD + FD + window (close to paper's estimator)
//   - Full: PSD + FD + window + bispectrum

package bt_config_pkg;

  // -----------------------------
  // Global numeric format
  // -----------------------------
  parameter int BT_WL        = 16;  // total bits
  parameter int BT_FRAC_BITS = 14;  // fractional bits (Q(1.FRAC))

  // -----------------------------
  // FFT size
  // -----------------------------
  parameter int BT_LOGN = 10;         // 2^10 = 1024
  parameter int BT_N    = (1 << BT_LOGN);

  // -----------------------------
  // Feature switches
  // -----------------------------
  // 1: enable fractional delay filter (6-tap FD of eq. (6)). 
  // 0: bypass FD (input goes straight to window/FFT).
  parameter bit ENABLE_FD = 1;

  // 1: enable windowing (Hann/Hamming/Blackman, etc.). 
  // 0: bypass window (rectangular).
  parameter bit ENABLE_WINDOW = 1;

  // 1: enable bispectrum B(f1,f2) = F(f1)F(f2)F*(f1+f2). 
  // 0: no bispectrum engine, PSD only.
  parameter bit ENABLE_BISPECTRUM = 0;

  // -----------------------------
  // Convenience "profiles"
  // -----------------------------
  // You can either set the above parameters manually when you
  // instantiate the top-level, or define compile-time profiles.
  //
  // Example:
  //   // Minimal PSD only:
  //   localparam bit P_MINIMAL_FD      = 0;
  //   localparam bit P_MINIMAL_WINDOW  = 0;
  //   localparam bit P_MINIMAL_BISPEC  = 0;
  //
  //   // Full paper-like:
  //   localparam bit P_FULL_FD         = 1;
  //   localparam bit P_FULL_WINDOW     = 1;
  //   localparam bit P_FULL_BISPEC     = 1;

endpackage
