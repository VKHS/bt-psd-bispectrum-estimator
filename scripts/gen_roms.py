# gen_roms.py
# Generate .mem files for:
#   - FFT twiddles: twiddle_cos.mem, twiddle_sin_neg.mem
#   - Window coefficients: window_hann.mem
#
# These are used by cru_shared.v and window_fir.v

import numpy as np

def quantize(x, frac_bits, data_w):
    scale = 1 << frac_bits
    v = np.round(x * scale).astype(int)
    # Saturate to signed range
    min_int = -(1 << (data_w - 1))
    max_int = (1 << (data_w - 1)) - 1
    v = np.clip(v, min_int, max_int)
    return v

def write_mem_hex(filename, vals, data_w):
    with open(filename, "w") as f:
        for v in vals:
            if v < 0:
                v = (1 << data_w) + v
            f.write(f"{v:0{data_w//4}X}\n")

def gen_twiddles(N=512, data_w=16, frac_w=14):
    k = np.arange(N)
    theta = 2 * np.pi * k / N
    cos_vals = np.cos(theta)
    sin_neg_vals = -np.sin(theta)  # for e^{-jθ}

    cos_q = quantize(cos_vals, frac_w, data_w)
    sin_q = quantize(sin_neg_vals, frac_w, data_w)

    write_mem_hex("twiddle_cos.mem", cos_q, data_w)
    write_mem_hex("twiddle_sin_neg.mem", sin_q, data_w)
    print("Generated twiddle_cos.mem and twiddle_sin_neg.mem")

def gen_window_hann(L=1024, data_w=16, frac_w=14):
    n = np.arange(L)
    w = 0.5 - 0.5 * np.cos(2 * np.pi * n / (L - 1))  # Hann window
    w_q = quantize(w, frac_w, data_w)
    write_mem_hex("window_hann.mem", w_q, data_w)
    print("Generated window_hann.mem")

if __name__ == "__main__":
    DATA_W = 16
    FRAC_W = 14
    FFT_N  = 512
    WIN_L  = 1024  # or 512, depending on your frame length

    gen_twiddles(N=FFT_N, data_w=DATA_W, frac_w=FRAC_W)
    gen_window_hann(L=WIN_L, data_w=DATA_W, frac_w=FRAC_W)
