#!/usr/bin/env bash
# Run the self-checking testbenches in the Vivado simulator (xsim).  Same sources, same golden vectors, same pass line
# ("TB_RESULT: PASS") as `make -C sim all` with Icarus Verilog.  The last four targets are the Model Composer pack:
# the record-and-replay wrapper (mc), the same wrapper as the include-free bundle that the Black Box imports (mcb), the
# test-signal generator (sg) and the in-fabric self-test (st).
#
#   source <Vivado install>/settings64.sh        # e.g. /tools/Xilinx/2025.2/Vivado/settings64.sh
#   bash sim/xsim_smoke.sh [cru bfly sq fft psd bsp mc mcb sg st]     # default: all ten
#
# NOT EXECUTED in the environment this repository was developed in (no Vivado there).  The flow follows the documented
# xvlog / xelab / xsim command lines (UG900); if a tool version rejects an option, the command that failed is printed.
# Everything is written to sim/work/xsim, which git ignores.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="$ROOT/sim/work/xsim"

for tool in xvlog xelab xsim; do
    command -v "$tool" >/dev/null 2>&1 || { echo "$tool not found: source <Vivado install>/settings64.sh first" >&2; exit 2; }
done

CRU="rtl/cru_rom.v rtl/cru_rotate.v"
BFLY="rtl/delay_line.v $CRU rtl/r8_butterfly.v"
FFT="rtl/ram_1r1w.v $BFLY rtl/fft512_mem.v"
PSD="$FFT rtl/mag_sq_odd.v rtl/psd_engine.v rtl/bt_psd_top.v"
BSP="$FFT rtl/bispectrum_engine.v rtl/bt_bispectrum_top.v"
MC="$PSD modelcomposer/hdl/bt_psd_mc.v"
MCB="modelcomposer/hdl/bt_psd_mc_bundle.v"
SG="modelcomposer/hdl/sig_gen.v modelcomposer/hdl/sig_rom.v"
ST="$MC $SG modelcomposer/hdl/bt_psd_selftest.v"

declare -A TOP=( [cru]=tb_cru_rotate [bfly]=tb_r8_butterfly [sq]=tb_mag_sq_odd [fft]=tb_fft512_mem [psd]=tb_psd_top [bsp]=tb_bispectrum_top
                 [mc]=tb_bt_psd_mc [mcb]=tb_bt_psd_mc [sg]=tb_sig_gen [st]=tb_bt_psd_selftest )
declare -A SRC=( [cru]="tb/tb_cru_rotate.v $CRU" [bfly]="tb/tb_r8_butterfly.v $BFLY" [sq]="tb/tb_mag_sq_odd.v rtl/mag_sq_odd.v"
                 [fft]="tb/tb_fft512_mem.v $FFT" [psd]="tb/tb_psd_top.v $PSD" [bsp]="tb/tb_bispectrum_top.v $BSP"
                 [mc]="tb/tb_bt_psd_mc.v $MC" [mcb]="tb/tb_bt_psd_mc.v $MCB" [sg]="tb/tb_sig_gen.v $SG" [st]="tb/tb_bt_psd_selftest.v $ST" )

TARGETS=("$@")
[ ${#TARGETS[@]} -gt 0 ] || TARGETS=(cru bfly sq fft psd bsp mc mcb sg st)

# The testbenches open "sim/vectors/*.mem" relative to the working directory: give the work directory a link to the vectors.
mkdir -p "$WORK/sim"
ln -sfn "$ROOT/sim/vectors" "$WORK/sim/vectors"
cd "$WORK"

status=0
for t in "${TARGETS[@]}"; do
    top="${TOP[$t]:?unknown target $t (cru bfly sq fft psd bsp mc mcb sg st)}"
    files=""
    for f in ${SRC[$t]}; do files="$files $ROOT/$f"; done
    echo "=== $t ($top)"
    # shellcheck disable=SC2086
    xvlog -i "$ROOT/rtl" -i "$ROOT/sim/vectors" $files > "$t.xvlog.log" 2>&1 || { echo "xvlog failed, see $WORK/$t.xvlog.log"; tail -n 15 "$t.xvlog.log"; status=1; continue; }
    xelab -debug off -timescale 1ns/1ps -s "${top}_sim" "work.$top" > "$t.xelab.log" 2>&1 || { echo "xelab failed, see $WORK/$t.xelab.log"; tail -n 15 "$t.xelab.log"; status=1; continue; }
    xsim "${top}_sim" -runall > "$t.log" 2>&1 || true
    grep -h "TB_RESULT" "$t.log" || { echo "TB_RESULT: FAIL $top (no result line, see $WORK/$t.log)"; status=1; continue; }
    grep -q "TB_RESULT: PASS" "$t.log" || status=1
done
exit $status
