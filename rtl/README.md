# RTL

Verilog-2001 sources of the reference implementation (design notes: `docs/ARCHITECTURE.md`).

**Status: simulated, not synthesised.** `make -C sim all` compiles the sources with Icarus Verilog (`-g2005 -Wall`) and
runs six self-checking testbenches against the golden vectors of the model; all six pass, bit for bit (CI repeats this
on every push). The sources also pass the structural lint (`scripts/lint_rtl.py`), and the control and addressing logic
has a line-by-line Python transliteration (`python/btpsd/rtlmodel.py`, `python/tests/test_rtlmodel.py`). The tested
configurations are `WL = 16` for the blocks and `WL = 12`, Hann window, `ROUND = 1` for the top levels; other word
lengths, windows and `ROUND = 0` need vectors generated for them. There is no synthesis, timing or power result.

| File | Content |
|---|---|
| `cru_rotate.v`, `cru_rom.v` | constant-coefficient rotator `y = x exp(-j q pi/512)`, 4 pipeline stages, generated ROM |
| `r8_butterfly.v` | radix-8 butterfly with output twiddle rotations, latency 10, one butterfly per cycle |
| `fft512_mem.v` | in-place 512-point radix-8 FFT on eight banked memories, load and unload ports, scale 1/512 |
| `mag_sq_odd.v` | magnitude squared, serial sum-of-odds squarer or multipliers |
| `psd_engine.v`, `bt_psd_top.v` | chunk merge (6-tap fractional delay), frequency-domain window, squarer, accumulation; top level with two FFT cores |
| `bispectrum_engine.v`, `bt_bispectrum_top.v` | window, rounded product, `P12 conj(F3)`, accumulation; top level with framing (hop 256) |
| `delay_line.v`, `ram_1r1w.v` | helpers: shift register, synchronous read-first RAM |
| `btpsd_params.vh` | latencies shared by the sources (`LAT_CRU`, `LAT_BFLY`) |
| `cru_rom.v`, `btpsd_consts.vh` | **generated** by `python/scripts/gen_rtl.py` (ROM, filter taps `FD_H0..5`, window taps `WIN_*`); do not edit |

Parameters worth knowing: `WL` (data word length; the testbenches use 16 for the blocks and 12 for the top levels, other values need vectors generated for that word length),
`ROUND` (1: round half up, 0: truncate), `ACCW` / `BACCW` (accumulator widths; worst-case bounds in
`docs/ARCHITECTURE.md`), `LOGB` (bispectrum bins per axis `2**LOGB`, at most 8), `SQ_SERIAL` (serial squarer).
The window defaults to Hann; the taps of the other windows are available as `WIN_RECT_*`, `WIN_HAMMING_*`,
`WIN_BLACKMAN_*` and are passed as parameters.

Regenerate the generated files and the vectors with `make -C sim vectors` (needs numpy and scipy).
