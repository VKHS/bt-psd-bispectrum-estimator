# Convenience targets. Needs python3 with numpy, scipy, matplotlib, pytest (see requirements.txt);
# `make sim` additionally needs Icarus Verilog (iverilog, vvp).
PY ?= python3

.PHONY: help test lint sim vectors experiments figures check

help:
	@echo "test         unit tests of the Python model and of the RTL transliteration"
	@echo "lint         structural check of the Verilog sources"
	@echo "sim          RTL testbenches with Icarus Verilog (needs the golden vectors in sim/vectors)"
	@echo "vectors      regenerate rtl/cru_rom.v, rtl/btpsd_consts.vh and sim/vectors/*"
	@echo "experiments  regenerate results/*.csv (about half a minute)"
	@echo "figures      regenerate docs/figures/*"
	@echo "check        test + lint"

test:
	$(PY) -m pytest -q

lint:
	$(PY) scripts/lint_rtl.py

check: test lint

sim:
	$(MAKE) -C sim -k all

vectors:
	$(MAKE) -C sim vectors

experiments:
	$(PY) python/scripts/run_experiments.py

figures:
	$(PY) python/scripts/make_diagrams_a.py
	$(PY) python/scripts/make_diagrams_b.py
	$(PY) python/scripts/make_diagrams_c.py
	$(PY) python/scripts/make_plots.py
