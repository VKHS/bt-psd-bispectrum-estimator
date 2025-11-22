### Overview

This repository collects the **software and hardware implementation** of a Blackman–Tukey power spectral density (PSD) and bispectral estimator designed for **non‑stationary biomedical signals** such as EEG, EOG and EMG. The design follows the architecture proposed in:

> A. V. Khalili Sadaghiani, B. Forouzandeh,
> *High‑performance power spectral/bispectral estimator for biomedical signal processing applications using novel memory‑based FFT processor*, Integration, the VLSI Journal, 2024. 

The core idea is to build a **resource‑efficient, low‑power estimator** that still gives high‑quality spectra when implemented with **short word lengths** on FPGA. To do this, the architecture combines several key techniques:

* A **memory‑based 1024‑point FFT** (here modeled at 512/1024 points) that reuses a single butterfly/rotational engine instead of a large parallel FFT, reducing area and power. 
* A **shared CORDIC‑II central rotational unit (CRU)** used not only for FFT twiddle rotations, but also by the fractional‑delay and window filters, so all expensive rotations live in one highly utilized block. 
* **Modified safe scaling** inside the FFT stages so that the FFT results are already scaled by (1/L). This avoids the explicit divisions and averaging usually needed at the end of a Blackman–Tukey estimator and makes fixed‑point implementation more robust. 
* A **6‑tap bidirectional fractional‑delay filter** to merge overlapping FFT “chunks” into the final Blackman–Tukey segments without having to compute a full extra FFT stage. 
* Optional **bispectrum computation** (B(f_1,f_2)=F(f_1)F(f_2)F^*(f_1+f_2)) for higher‑order analysis of nonlinearities and phase couplings in biomedical signals. 

The repository is organised into three layers:

1. **Python models** – floating‑point and fixed‑point reference implementations of the BT PSD and bispectrum, plus a structural model that mirrors the RTL.
2. **Core HDL modules** – synthesizable Verilog/SystemVerilog blocks: memory‑based FFT, CRU, fractional‑delay filter, window filter, magnitude‑squared and accumulation, bispectrum engine, and ROM generators.
3. **Top‑level integration** – configuration package and a configurable top module that can be built as:

   * a **minimal PSD** engine (FFT + |X|²),
   * a **full BT PSD** estimator with FD + windowing, or
   * a **full PSD + bispectrum** system.

Together these pieces let you go from **paper → Python → FPGA**: you can validate algorithms against the paper’s results, tune fixed‑point formats and then map the same architecture onto an FPGA for real‑time biomedical signal processing experiments.



---

## Repository structure

```text
.
├── 1. Three main code bases/        # Python / reference code
│   ├── bt_reference.py              # End-to-end BT PSD/bispectrum reference
│   ├── bt_structural.py             # Structural model mirroring HDL blocks
│   ├── fft_cordic_fixed.py          # Fixed-point CORDIC FFT reference
│   └── fixed_point_model.py         # Helpers for quantisation and word-length sweeps
│
├── 2. Core HDL modules/             # Reusable Verilog/SystemVerilog building blocks
│   ├── complex_fixed.v              # Complex add/sub/mul in fixed point
│   ├── cru_shared.v                 # Shared CORDIC-II rotational unit (CRU)
│   ├── dual_port_ram(updated).sv    # Dual-port RAM used by memory-based FFT
│   ├── fft_radix2_mem.v             # Memory-based radix-2 FFT (reference variant)
│   ├── fft512_mem_core.sv           # 512-pt FFT core with safe-scaling
│   ├── twiddle_rom.sv               # Twiddle ROM (cos/sin) for FFT
│   ├── fd6_fir.v                    # 6-tap bidirectional fractional-delay filter (FD)
│   ├── window_fir.v                 # Folded window filter (Hann/Hamming/Blackman)
│   ├── mag_sq_odd.v                 # Sum-of-odds squarer (multiplierless |X|²)
│   ├── psd_accumulator.v            # Optional PSD averaging / save-add unit
│   ├── bispectrum_engine.v          # Triple product core B(f1,f2)=F(f1)F(f2)F*(f1+f2)
│   ├── gen_twiddle_rom.py           # Python script to generate FFT twiddle .mem files
│   └── gen_roms.py                  # Any additional ROM generators
│
└── 3. Additional necessary/         # Config + top-level integration
    ├── bt_config_pkg.sv             # Global parameters + feature switches
    └── bt_estimator_top.sv          # Full BT estimator top-level (FFT + FD + window + PSD + optional bispectrum)
