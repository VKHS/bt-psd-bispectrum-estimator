# gen_twiddle_rom.py
#
# Generate twiddle ROM init files:
#   twiddle_re_512_q15.mem
#   twiddle_im_512_q15.mem
#
# For W_N^k = exp(-j*2*pi*k/N), k = 0..N/2-1, in Q1.15 format.

import numpy as np

N        = 512
DATA_W   = 16
FRAC_W   = 15
SCALE    = 1 << FRAC_W
K_MAX    = N // 2  # 256 entries

def to_q15(x):
    """Quantize float -> signed Q1.15 integer."""
    v = np.round(x * SCALE).astype(int)
    min_int = -(1 << (DATA_W - 1))
    max_int = (1 << (DATA_W - 1)) - 1
    v = np.clip(v, min_int, max_int)
    return v

def to_hex_twos_complement(v):
    if v < 0:
        v = (1 << DATA_W) + v
    return f"{v:0{DATA_W//4}X}"

k = np.arange(K_MAX)
theta = -2.0 * np.pi * k / N  # negative sign for FFT

re = np.cos(theta)
im = np.sin(theta)

re_q = to_q15(re)
im_q = to_q15(im)

with open("twiddle_re_512_q15.mem", "w") as f_re, \
     open("twiddle_im_512_q15.mem", "w") as f_im:
    for r, i in zip(re_q, im_q):
        f_re.write(to_hex_twos_complement(r) + "\n")
        f_im.write(to_hex_twos_complement(i) + "\n")

print("Generated twiddle_re_512_q15.mem and twiddle_im_512_q15.mem")
