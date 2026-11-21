# Convenience targets. Needs python3 with numpy, scipy, matplotlib, pillow, pytest (see requirements.txt);
# `make sim` additionally needs Icarus Verilog (iverilog, vvp).
PY ?= python3

.PHONY: help test lint smoke sim matrix vectors experiments figures check

help:
	@echo "test         unit tests of the Python model and of the RTL transliteration"
	@echo "lint         structural check of the Verilog sources"
	@echo "smoke        quick check of the whole repository: files, generated sources, vectors, models, Model Composer pack, lint, figures; the testbenches too if iverilog is installed"
	@echo "sim          RTL testbenches with Icarus Verilog (needs the golden vectors in sim/vectors)"
	@echo "matrix       RTL against the model for 25 word-length, window, rounding and squarer variants (Icarus Verilog)"
	@echo "vectors      regenerate rtl/cru_rom.v, rtl/btpsd_consts.vh and sim/vectors/*"
	@echo "experiments  regenerate results/*.csv (about half a minute)"
	@echo "figures      regenerate docs/figures/* (31 sheets; fig19, fig20, fig23 and fig31 read results/rtl_*.csv)"
	@echo "check        test + lint"

test:
	$(PY) -m pytest -q

lint:
	$(PY) scripts/lint_rtl.py

smoke:
	$(PY) python/scripts/smoke_test.py

check: test lint

sim:
	$(MAKE) -C sim -k all

matrix:
	$(MAKE) -C sim matrix

vectors:
	$(MAKE) -C sim vectors

experiments:
	$(PY) python/scripts/run_experiments.py

figures:
	$(PY) python/scripts/make_diagrams_a.py
	$(PY) python/scripts/make_diagrams_b.py
	$(PY) python/scripts/make_diagrams_c.py
	$(PY) python/scripts/make_plots.py
	$(PY) python/scripts/make_rtl_figures.py
	$(PY) python/scripts/make_sim_figures.py
	$(PY) python/scripts/make_signal_figures.py
	$(PY) python/scripts/make_fft_figures.py
	$(PY) python/scripts/make_spectral_figures.py
	$(PY) python/scripts/make_mc_figures.py
