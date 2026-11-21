# Model Composer pack

The PSD estimator as a Black Box for AMD Vitis Model Composer (formerly System Generator for DSP), a Simulink test bench with
signal-generator sources around it, and a test-signal generator that runs inside the FPGA fabric without MATLAB. The full guide,
with the block-by-block bench, the expected results and the clock-cycle counts, is [`docs/MODEL_COMPOSER.md`](../docs/MODEL_COMPOSER.md).

| File | Role |
|---|---|
| `hdl/bt_psd_mc.v` | the Black Box wrapper: records the samples, replays them into `rtl/bt_psd_top.v`, raises `done` |
| `hdl/bt_psd_mc_bundle.v` | the wrapper and all the modules it needs in one include-free file: the file the Black Box compiles (generated) |
| `hdl/sig_gen.v`, `hdl/sig_rom.v` | in-fabric test-signal generator: phase accumulator, sine table, LFSR noise (`sig_rom.v` generated) |
| `hdl/bt_psd_selftest.v` | generator and wrapper together: the estimator on a generated tone |
| `matlab/bt_psd_mc_config.m` | Black Box configuration function: ports, types, clock and clock enable, HDL file (generated) |
| `matlab/bt_mc_params.m` | word length, chunks, segments, averaging shift, widths and power scale of the build (generated) |
| `matlab/bt_mc_build.m` | builds the Simulink bench: Sine Wave and Random Number, Signal Generator, Chirp Signal or From Workspace into Gateway In |
| `matlab/bt_mc_run.m` | simulates the bench, saves the logs, runs the checker |

The generated files come from `python/scripts/gen_mc.py`; `python/tests/test_gen_mc.py` fails when they are out of date. The
folders `stimulus/` and `work/` receive what `export_mc_stimulus.py` and the MATLAB scripts write and are ignored by git.

## Two routes

**Simulink.** Start MATLAB through Model Composer, then from the repository root:

```matlab
>> addpath('modelcomposer/matlab')
>> mdl = bt_mc_build('sine');      % 'signal', 'chirp' or ('file', 'File', <stimulus .mat>)
>> res = bt_mc_run(mdl);           % simulates, then checks the spectrum word for word against the Python model
```

**Vivado simulator or Icarus Verilog, no MATLAB.** The self-test drives the estimator from `sig_gen.v` and compares the spectrum
with the model:

```bash
make -C sim sg st mc mcb                     # Icarus Verilog
bash sim/xsim_smoke.sh sg st mc mcb          # Vivado simulator, after sourcing settings64.sh of the installation
```

## Status

The Verilog is simulated with Icarus Verilog (the CI workflow runs it) and agrees with the bit-exact Python model; the clock
cycles of the simulation equal the model's prediction on the testbench records. The MATLAB scripts and `sim/xsim_smoke.sh` are
**not executed** where this repository was written, because there is no MATLAB, Simulink, Model Composer or Vivado there; they
are checked statically against the Verilog and the Python tools. The guide lists what is established and what is not. Nothing in
this pack is synthesised, and no area, frequency or power figure is claimed.
