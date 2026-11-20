# Original files

The files of the repository as it was before the expansion, unchanged and in their original folders:

| Folder | Content |
|---|---|
| `1. Three main code bases/` | `bt_reference.py`, `bt_structural.py`, `fft_cordic_fixed.py`, `fixed_point_model.py` |
| `2. Core HDL modules/` | rotator, FIR, FFT, squarer, accumulator and bispectrum modules, ROM generators |
| `3. Additional necessary/` | `bt_config_pkg.sv`, `bt_estimator_top.sv` |

Eight of the files carry stray `:contentReference[oaicite:N]{index=N}` markers at the end of comment and docstring lines (an artefact of a text export, not part of the code); they are left as they were.

They are kept for reference only. The maintained code is in `python/` (models), `rtl/` (RTL), `tb/` and `sim/`
(testbenches and golden vectors); none of the tests, the lint or the simulation build reads this folder.
