# Notes on the paper: findings from building the model and candidate errata

Paper: A. Khalili Sadaghiani, B. Forouzandeh, *High-performance power spectral/bispectral estimator for biomedical
signal processing applications using novel memory-based FFT processor*, Integration, the VLSI Journal 99 (2024) 102241.

Values marked *paper* are quoted from the paper (`data/paper_*.csv`). Values marked *model* were measured in this
repository with `python/scripts/run_experiments.py` on synthetic signals (`results/*.csv`); the paper's EEG, EOG and
EMG records are not part of the repository, so *model* and *paper* values are never merged into one number.

## 1. Findings that matter for anyone using the numbers

| # | Finding | Evidence | Consequence |
|---|---|---|---|
| 1 | The estimator computes **Welch's method** (averaged windowed periodograms), not a Blackman-Tukey estimate | section 2.1 | the reference in the code is the Welch average; the name `bt` is kept for the repository only |
| 2 | The 6-tap fractional-delay merge limits the odd-bin accuracy to about **11.7 dB**; Table 3 reports 50.08 dB for six taps and the model does not reproduce it | section 2.2, `results/fd_merge_snr_vs_taps.csv` | the PSD SQNR of every word length saturates at a merge limit (23.5 to 44 dB, depending on the signal) |
| 3 | The SQNR of Table 10 is a single-point value at (f1, f2) = (10, 8); the model reports plane-wide values | section 2.3 | no head-to-head comparison of Table 10 with the model |
| 4 | The sum-of-odds squarer needs `max(|re|, |im|)` cycles; at WL = 16 one segment needs 25 k to 56 k squarer cycles | section 2.4, `results/squarer_cycle_model.csv` | a serial squarer meets the throughput of Table 9 at WL = 10 but not at WL = 16 |
| 5 | Twelve candidate errata and omissions in equations and tables | section 3 | listed for the authors' own check against the source files |

## 2. Details

### 2.1 Welch, not Blackman-Tukey

The merge produces the 1024-point spectrum `X(k)` of a data segment. The window is applied to the *complex* spectrum
(Eq. 13, `Xw(k) = sum_d W(d) X(k-d)`), where `W(d)` are the DFT coefficients of a data window `w[n]`
(exact for the raised-cosine windows, see `btpsd.windows.spectral_kernel`). Hence `Xw = DFT(w . x)`, and the average of
`|Xw|^2` over overlapping segments is the Welch estimate:

```
P_Welch(k) = (1/K) sum_s |DFT(w . x_s)(k)|^2
```

A Blackman-Tukey estimate multiplies the *autocorrelation estimate* by a lag window and transforms it, which equals
a convolution of the periodogram `|X|^2` (a real, non-negative sequence) with the transform of the lag window. The
repository therefore uses `btpsd.reference.welch_twosided(scaling="paper")` as the exact reference. The bispectrum is
the direct (segment-averaged triple product) estimator, with `bispectrum_indirect` (third-order cumulant) available
for comparison.

### 2.2 The fractional-delay merge

For a segment made of chunks `a[n]` and `b[n]` (512 samples each) the odd bins are, exactly,

```
X(2m+1) = sum_n (a[n] - b[n]) e^{-j pi n/512} e^{-j 2 pi m n/512} = sum_d g(d) D(m + d),    D = DFT_512(a - b)
|g(d)| = 1 / (pi |d + 1/2|)    ->    0.6366, 0.2122, 0.1273, 0.0909, ...      sum over all d of |g|^2 = 1
```

Truncating to `2K` taps keeps `E_K = sum_{kept} h^2` of the kernel energy and discards the rest, so for noise-like or
line spectra the error-to-signal ratio of the odd bins is `1 - E_K` (`E_K = sum_{kept} |g|^2`):

| taps | kept energy | bound `-10 log10(1 - E_K)` | model (white / pink / brown / EEG-like / AR(2)) | paper, Table 3 |
|---|---|---|---|---|
| 4 | 0.9006 | 10.03 dB | 10.0 / 10.1 / 9.7 / 10.1 / 10.1 dB | 37.86 dB |
| 6 | 0.9331 | 11.74 dB | 11.7 / 11.9 / 11.4 / 11.9 / 11.8 dB | 50.08 dB |
| 8 | 0.9496 | 12.98 dB | 12.9 / 13.2 / 12.6 / 13.2 / 13.1 dB | 51.46 dB |
| 10 | 0.9596 | 13.94 dB | 13.8 / 14.2 / 13.6 / 14.2 / 14.0 dB | 51.97 dB |

The six taps of Eq. 6 (third tap 0.125 instead of 0.1273) keep 0.9318, i.e. 11.66 dB. The model agrees with the bound
to within 0.4 dB for five different signals; Table 3 is 28 to 38 dB above it for even tap counts, and its odd tap
counts (8.14, 15.21, 20.67, 22.13 dB) lie 23 to 35 dB below their even neighbours, a gap the model does not show (its
SNR rises smoothly from 8.4 dB at three taps to 13.8 dB at ten for white noise). The paper does not define the SNR of
Table 3 or name its test signal. If it is computed on something other than the complex odd-bin spectrum (only on
smooth parts of the spectrum, for example), the definition is what is needed to reproduce it; with a definition and
a record, `run_experiments.py` can be pointed at them.

The merge limit propagates to the PSD: against the exact Welch average the model saturates at 23.55 dB (white),
28.84 dB (pink), 44.01 dB (EEG-like) and 23.34 dB (AR(2)) as the word length grows.

### 2.3 What the SQNR figures measure

* *Paper:* Figs. 18-20 compare the FPGA output with "64-bit MATLAB" periodogram and Blackman-Tukey PSDs, and the SQNR
  of Tables 8 and 10 is computed against IEEE-754 64-bit results. The formula is not given.
* *Model:* `btpsd.metrics.sqnr_db` = `10 log10(sum P_ref^2 / sum (P_ref - P_hw)^2)` on linear power at the analytic scale, reported
  against two references: the exact Welch average (includes the merge error) and the floating-point six-tap merge
  (quantisation error only).
* *Table 10 (paper):* SQNR "for the case of f1 = 10, and f2 = 8", three cases (tau1 = tau2 = 0.85; tau1 = tau2 = 1.3;
  tau1 = 1, tau2 = 1.75, rectangular window). The role of time lags tau in a frequency-domain estimator is not
  explained. The values (8.6 to 21.4 dB for WL = 10 to 18) refer to one bispectrum value.
* *Model, bispectrum, magnitude SQNR over the whole 128 x 128 plane (coupled quadratic-phase-coupling signal):*
  22.7, 34.6, 43.1, 52.2, 59.6 dB at WL = 10, 12, 14, 16, 18 (`results/bispectrum_sqnr_vs_wordlength.csv`).
* *Model, PSD SQNR against the exact Welch average, WL = 10 / 12 / 14 / 16 / 18 / 20:*

| signal | 10 | 12 | 14 | 16 | 18 | 20 | merge limit |
|---|---|---|---|---|---|---|---|
| white | 15.4 | 23.9 | 23.6 | 23.5 | 23.5 | 23.5 | 23.55 |
| pink | 28.4 | 31.9 | 29.5 | 29.0 | 28.9 | 28.9 | 28.84 |
| EEG-like | 28.7 | 40.0 | 44.2 | 44.2 | 44.0 | 44.0 | 44.01 |
| AR(2) narrow-band | 21.6 | 23.4 | 23.3 | 23.3 | 23.3 | 23.3 | 23.34 |

From WL = 16 the SQNR is within 0.2 dB of the merge limit. At WL = 12 and 14 the white and pink values lie above it (by up
to 3 dB for pink at WL = 12); this was not investigated further. Table 8 (*paper*, WL = 10, proposed hardware) lists 11.8 to 20.6 dB for the nine records, the
same range as the model at WL = 10 (15.4 to 28.7 dB), but the records differ and no closer comparison is made.

### 2.4 The sum-of-odds squarer and throughput

`n^2 = 1 + 3 + ... + (2|n| - 1)` takes `|n|` cycles, so the cycle count follows the sample code (two squarers in
parallel: `max(|re|, |im|)`). Model results for 513 bins per segment (`results/squarer_cycle_model.csv`):

| signal | WL | mean cycles / segment | worst bin | squarer pairs needed for R = 0.88 |
|---|---|---|---|---|
| white | 10 | 919 | 7 | 1.58 |
| white | 16 | 55 723 | 369 | 95.8 |
| pink | 10 | 702 | 71 | 1.21 |
| pink | 16 | 40 113 | 4 462 | 68.9 |
| EEG-like | 10 | 606 | 44 | 1.04 |
| EEG-like | 16 | 31 332 | 2 780 | 53.9 |
| AR(2) | 10 | 576 | 23 | 0.99 |
| AR(2) | 16 | 25 633 | 1 466 | 44.1 |

"Pairs needed" is the mean cycle count divided by the 582 cycles that a chunk of 512 samples may take at the
throughput ratio R = 0.88 samples per cycle of Table 9. At WL = 10 one or two serial squarer pairs suffice; at WL = 16
tens of pairs are needed or the squarer must be a multiplier. The text of the paper shifts the value "to the left by half
the size of WL, in order to make it an integer number" before counting. This is read here as aligning the binary
point, so the count equals the integer sample code. If it instead means counting only the upper `WL/2` bits, the count
drops to at most `2^(WL/2)` cycles but the squared value loses precision; the paper does not say which.

## 3. Candidate errata and omissions

Each item was checked on the printed page. They are candidates for the authors to confirm against the source files.

| # | Where | As printed | Expected |
|---|---|---|---|
| 1 | Eq. (2) | outer sum `F = 0..511` with samples `x(m + 8n + 64l + 2F)` | the largest index is 7 + 56 + 448 + 1022 = 1533, beyond a 1024-sample segment; the three inner sums already form a 512-point DFT, so the outer index takes two values (for example `x(2u + F)`, `F = 0, 1`, `u = m + 8n + 64l`) |
| 2 | Table 5, row `tg(w_i) = 2^n` | 1, 0.25, 0.125, 0.0625, 0.03125, 0.03125, 0.01562, 0.00781 | `2^-n`: 1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125 (entries n = 1..4 are shifted by one; the angle row is right for n = 0..4) |
| 3 | Table 5, angle row n = 5, 6, 7 | 1.787, 0.893, 0.446 degrees | 1.790, 0.895, 0.448 degrees |
| 4 | text before Eq. (18) | `x1 = x0 - y0 + tg(w0)`, `y1 = y0 - x0 + tg(w0)` | `x1 = x0 - y0 tg(w0)`, `y1 = y0 + x0 tg(w0)` (as in the step to `x2, y2`) |
| 5 | Eq. (17) | `[(x0 - y0 tg(w0)) + (y0 + x0 tg(w0))]` | the imaginary unit is missing: `(x0 - y0 tg(w0)) + j (y0 + x0 tg(w0))` |
| 6 | Eq. (19) | `(x2 + j y2) e^{j w1} e^{j w2} ...` | `w1` has been applied; the remaining product starts at `e^{j w2}` |
| 7 | Eqs. (20)-(23) | (21) and (23) define `x_{i+1}` a second time; (20), (21) use `tg(w0)` | `y_{i+1}` in (21) and (23); `tg(w_i)` in (20), (21) |
| 8 | Eq. (30) | `|B| = |F(f1)| + |F(f2)| + |F(f1+f2)|` | `B = F(f1) F(f2) F*(f1+f2)` gives the product of the three magnitudes (the phase relation (29) is the one with sums) |
| 9 | Eq. (32) | third-order cumulant with products of `C2y` terms and a third lag `tau3` | `C2y` and `tau3` are not defined for a two-lag cumulant; the `C2y C2y` pairings are the terms of the fourth-order cumulant. The repository uses the segment-averaged third-order moment of mean-removed data (`third_order_cumulant`) |
| 10 | Table 6 | 128 rows with a label `k pi/512` and a coefficient pair `(a, b)` | in 27 of 128 rows the angle of `(a, b)` differs from the label by more than the stated 0.01 rad (1.63 steps of pi/512), up to 33 steps (row 8: label pi/64 with (20, 5)); `data/table6_printed.csv` lists the rows and `paperdata.table6_printed()` flags them. The ROM of this repository does not use Table 6 |
| 11 | Table 2, row R113 | first entry 340, later entry 1361 | the neighbours follow `8 i - 504` (R112: 392, R114: 408) and +32 steps (1505, 1633), which gives 400 and 1601 |
| 12 | Table 8 | columns "Unrolled C." and "Single C." | no values under either heading; the caption lists the window order Rectangular, Hamming, Hann, Blackman, which matches the row order |
| 13 | Section 2 (page 3) against Fig. 3 | the text calls the barrel shifter fifteen-bit three times; Fig. 3 labels it "8-bit Barrel Shifter (For choosing Ms)" inside a 15-bit main counter made of a 7-bit pass counter and an 8-bit butterfly counter | one width; the counter split 7 + 8 = 15 is consistent, the shifter width is not |
| 14 | Section 2 (page 3) | "Eight memory-banks in eight clock cycles receive sixty-eight samples and send sixty-eight samples. This is also shown in Tables 2 and 3." | 8 banks x 8 cycles = 64 words; the register-sequence tables are Tables 1 and 2 (Table 3 is the SNR of the fractional-delay taps) |

Item 3 changes nothing in the results; items 1, 4-7 matter if the equations are used as implementation references.

## 4. What the paper leaves open (and what this repository chose)

| Open point | Choice here |
|---|---|
| rounding of every operation | round half up, saturation only where overflow is possible |
| fixed-point format of the six taps and window taps | 12 and 10 fraction bits (`FD_FRAC`, `WIN_FRAC`) |
| CRU coefficients | generated ROM, 129 entries, 0.01 rad tolerance (Table 6 is not self-consistent, item 10) |
| accumulator widths | `2 WL + 8` (PSD), `2 WL + 10` (bispectrum) with the worst-case bounds in `ARCHITECTURE.md` |
| SNR definition of Table 3 and the signal | not defined; see section 2.2 |
| meaning of tau in Table 10 | not defined; see section 2.3 |
| squarer count (shifted code or full code) | full integer code, see section 2.4 |
