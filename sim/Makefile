# RTL simulation with Icarus Verilog.  Run `make -C sim all` (or one of: cru bfly sq fft psd bsp mc mcb sg st) from anywhere.
# `make -C sim matrix` runs the PSD and bispectrum top levels for 25 further configurations (needs numpy and scipy).
# The testbenches read sim/vectors/*.mem relative to the repository root, so every run starts there.
# Regenerate the vectors and the generated RTL constants with `make -C sim vectors` (needs numpy and scipy).
ROOT   := $(abspath ..)
IV     ?= iverilog
VVP    ?= vvp
IVFLAGS = -g2005 -Wall -I$(ROOT)/rtl -I$(ROOT)/sim/vectors
WORK   := $(ROOT)/sim/work

CRU  := rtl/cru_rom.v rtl/cru_rotate.v
BFLY := rtl/delay_line.v $(CRU) rtl/r8_butterfly.v
FFT  := rtl/ram_1r1w.v $(BFLY) rtl/fft512_mem.v
PSD  := $(FFT) rtl/mag_sq_odd.v rtl/psd_engine.v rtl/bt_psd_top.v
BSP  := $(FFT) rtl/bispectrum_engine.v rtl/bt_bispectrum_top.v
MC   := $(PSD) modelcomposer/hdl/bt_psd_mc.v
MCB  := modelcomposer/hdl/bt_psd_mc_bundle.v
SG   := modelcomposer/hdl/sig_gen.v modelcomposer/hdl/sig_rom.v
ST   := $(MC) $(SG) modelcomposer/hdl/bt_psd_selftest.v

.PHONY: all cru bfly sq fft psd bsp mc mcb sg st matrix vectors clean
all: cru bfly sq fft psd bsp mc mcb sg st
	@echo "--- summary ---"; grep -h "TB_RESULT" $(foreach t,cru bfly sq fft psd bsp mc mcb sg st,$(WORK)/$(t).log)

define run
	@mkdir -p $(WORK)
	cd $(ROOT) && $(IV) $(IVFLAGS) -o sim/work/$(1).vvp $(2) && $(VVP) sim/work/$(1).vvp | tee sim/work/$(1).log
	@grep -q "TB_RESULT: PASS" $(WORK)/$(1).log
endef

cru:  ; $(call run,cru,tb/tb_cru_rotate.v $(CRU))
bfly: ; $(call run,bfly,tb/tb_r8_butterfly.v $(BFLY))
sq:   ; $(call run,sq,tb/tb_mag_sq_odd.v rtl/mag_sq_odd.v)
fft:  ; $(call run,fft,tb/tb_fft512_mem.v $(FFT))
psd:  ; $(call run,psd,tb/tb_psd_top.v $(PSD))
bsp:  ; $(call run,bsp,tb/tb_bispectrum_top.v $(BSP))
mc:   ; $(call run,mc,tb/tb_bt_psd_mc.v $(MC))
mcb:  ; $(call run,mcb,tb/tb_bt_psd_mc.v $(MCB))
sg:   ; $(call run,sg,tb/tb_sig_gen.v $(SG))
st:   ; $(call run,st,tb/tb_bt_psd_selftest.v $(ST))

matrix:
	cd $(ROOT) && PYTHONPATH=python python3 python/scripts/run_rtl_matrix.py --iverilog $(IV) --vvp $(VVP)

vectors:
	cd $(ROOT) && PYTHONPATH=python python3 python/scripts/gen_rtl.py && PYTHONPATH=python python3 python/scripts/gen_vectors.py

clean:
	rm -rf $(WORK)
