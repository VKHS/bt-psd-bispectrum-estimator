# Figures of the paper

These 32 images are the figures of

> AbdolVahab Khalili Sadaghiani and Behjat Forouzandeh, "High-performance power spectral/bispectral estimator for biomedical signal
> processing applications using novel memory-based FFT processor", *Integration, the VLSI Journal* 99 (2024) 102241.

They are **paper-reported**: they were produced by the authors for the article and were **not produced by this
repository**. Nothing in `docs/figures/` (the sheets drawn by `python/scripts/make_*.py`) is copied from them.
The paper's plots of recorded EEG, EOG and EMG signals (Figs. 18-32) cannot be regenerated here because the
recordings are not part of the repository; the model's counterparts use synthetic signals and say so in
their titles.

Rights: the repository carries no licence file, and nothing in it grants rights to these images. The published article
carries an Elsevier B.V. copyright notice (2024) that reserves all rights, including text and data mining; the re-use of
its figures is governed by the publisher's author-rights terms. Check them before keeping or redistributing the images in a
public repository.

File names are `paperfigNN_<slug>`, `NN` being the figure number in the article. "Page" is the page of the
published PDF that carries the caption.

| Fig. | File | Page | Content | Counterpart in this repository |
|---:|---|---:|---|---|
| 1 | `paperfig01_spectral_estimation_taxonomy.png` | 2 | Spectral estimation of stochastic and deterministic signals | `README.md`, `docs/PAPER_NOTES.md` |
| 2 | `paperfig02_r8mdc_fft_pipeline.png` | 3 | 512-point R8MDC FFT pipelined processor (the baseline the memory-based design replaces) | `docs/figures/fig03_r8_butterfly_network` |
| 3 | `paperfig03_proposed_estimator_block_diagram.jpg` | 6 | Block diagram of the proposed estimator | `docs/figures/fig01_system_overview`, `fig02_memory_fft_architecture`, `fig15_rtl_hierarchy` (the Verilog hierarchy) |
| 4 | `paperfig04_chunk_segments_and_combination.png` | 7 | (a) Blackman-Tuckey segments from previous and next chunks, (b) the combination operation | `docs/figures/fig04_segment_merge` |
| 5 | `paperfig05_bidirectional_estimator.png` | 7 | Sub-sample bidirectional estimator | `docs/figures/fig05_fd6_filter` |
| 6 | `paperfig06_fd3_impulse_group_delay.png` | 7 | Impulse and group-delay response, 3-tap fractional delay | `docs/figures/fig11_fd_kernel_gallery` |
| 7 | `paperfig07_fd4_impulse_group_delay.png` | 7 | 4-tap fractional delay | `fig11_fd_kernel_gallery` |
| 8 | `paperfig08_fd5_impulse_group_delay.png` | 8 | 5-tap fractional delay | `fig11_fd_kernel_gallery` |
| 9 | `paperfig09_fd6_impulse_group_delay.jpg` | 8 | 6-tap fractional delay (the tap set used by the design) | `fig05_fd6_filter`, `fig11_fd_kernel_gallery` |
| 10 | `paperfig10_fd7_impulse_group_delay.png` | 8 | 7-tap fractional delay | `fig11_fd_kernel_gallery` |
| 11 | `paperfig11_fd8_impulse_group_delay.png` | 8 | 8-tap fractional delay | `fig11_fd_kernel_gallery` |
| 12 | `paperfig12_window_rectangular.png` | 9 | Frequency and impulse response, rectangular window | `fig21_window_responses`, `btpsd.windows`, `fig07_window_squarer_accumulator` |
| 13 | `paperfig13_window_hann.png` | 9 | Hann window | `fig21_window_responses` |
| 14 | `paperfig14_window_hamming.png` | 9 | Hamming window | `fig21_window_responses` |
| 15 | `paperfig15_window_blackman.png` | 10 | Blackman window | `fig21_window_responses` |
| 16 | `paperfig16_rotation_angle_choices_mcr_scr.jpg` | 10 | Two possible choices of rotation angles (MCR, SCR) | `fig06_cru_rotator` |
| 17 | `paperfig17_cru_block_diagram.png` | 10 | Block diagram of the central rotational unit | `fig06_cru_rotator`, `fig18_pipeline_timing` (a), `rtl/cru_rotate.v` |
| 18 | `paperfig18_psd_eeg1.jpg` | 16 | EEG-1 signal and its estimated Blackman-Tuckey PSD | `fig22_psd_examples` (synthetic signals), `fig12_wordlength_study` |
| 19 | `paperfig19_psd_eog1.jpg` | 16 | EOG-1 signal and its PSD | `fig22_psd_examples`, `fig12_wordlength_study` |
| 20 | `paperfig20_psd_emg1.jpg` | 16 | EMG-1 signal and its PSD | `fig22_psd_examples`, `fig12_wordlength_study` |
| 21 | `paperfig21_sqnr_eeg1.jpg` | 16 | SQNR of the proposed and competing methods, EEG-1 | `fig12_wordlength_study` |
| 22 | `paperfig22_sqnr_eog1.jpg` | 17 | SQNR, EOG-1 | `fig12_wordlength_study` |
| 23 | `paperfig23_sqnr_emg1.jpg` | 17 | SQNR, EMG-1 | `fig12_wordlength_study` |
| 24 | `paperfig24_bispectrum_eeg1_rectangular.jpg` | 17 | Bispectrum of EEG-1, rectangular window | `fig13_bispectrum_atlas`, `fig14_bispectrum_views` |
| 25 | `paperfig25_bispectrum_eeg1_hann.jpg` | 18 | EEG-1, Hann window | `fig13_bispectrum_atlas` |
| 26 | `paperfig26_bispectrum_eeg1_hamming.jpg` | 18 | EEG-1, Hamming window | `fig13_bispectrum_atlas` |
| 27 | `paperfig27_bispectrum_eog1_rectangular.jpg` | 18 | EOG-1, rectangular window | `fig13_bispectrum_atlas` |
| 28 | `paperfig28_bispectrum_eog1_hann.jpg` | 18 | EOG-1, Hann window | `fig13_bispectrum_atlas` |
| 29 | `paperfig29_bispectrum_eog1_hamming.jpg` | 19 | EOG-1, Hamming window | `fig13_bispectrum_atlas` |
| 30 | `paperfig30_bispectrum_emg1_rectangular.jpg` | 19 | EMG-1, rectangular window | `fig13_bispectrum_atlas` |
| 31 | `paperfig31_bispectrum_emg1_hann.jpg` | 19 | EMG-1, Hann window | `fig13_bispectrum_atlas` |
| 32 | `paperfig32_bispectrum_emg1_hamming.jpg` | 19 | EMG-1, Hamming window | `fig13_bispectrum_atlas` |

## How to read the two sets side by side

* Paper values drawn in this repository's plots are hollow black diamonds and are labelled "paper" or "Table 3";
  values measured on the model use filled markers. The two are never merged into one number.
* Figs. 24-32 show cumulant, magnitude and phase panels of the bispectrum of recorded signals for three windows.
  `fig13_bispectrum_atlas` repeats the *window* comparison on three synthetic signal classes (coupled EEG-like,
  quadratic phase coupling, pink noise) at a fixed word length, with zoomed views of the coupling peaks;
  `fig14_bispectrum_views` shows the cumulant, magnitude and phase views of one signal.
* The paper does not define the SNR behind its Table 3 and Figs. 21-23. The model's SQNR is defined in
  `btpsd.metrics.sqnr_db` and in the sheet subtitles; see `docs/PAPER_NOTES.md` for the comparison.
* Figs. 13-15 show *symmetric* 128-point windows: the heights of the printed stems agree with the symmetric windows and not
  with the periodic ones (`data/paper_figs13_15_stems_read.csv`, measured on the pixels by
  `python/scripts/measure_paper_window_stems.py`; `docs/PAPER_NOTES.md`, item 15). `fig21_window_responses` uses the periodic
  form of the datapath, for which the 1, 3 and 5-tap kernels are exact.
* `fig22_psd_examples` shows two of the paper's three PSD curves (the float periodogram and the fixed-point estimator at
  WL = 10); the third, an ARMA(18,20) estimator, is not part of this repository.
