# Out-of-context synthesis for a first resource and timing estimate (no I/O constraints, no place and route).
#
#   vivado -mode batch -source fpga/synth_ooc.tcl -tclargs <top> [part] [clock_ns] [WL]
#
#   top       bt_psd_top | bt_bispectrum_top | fft512_mem
#   part      default xc7a100tcsg324-1 (an Artix-7 device; the paper names only "Artix-7")
#   clock_ns  default 10.0
#   WL        data word length, default 12
#
# Run from the repository root. Reports and the checkpoint go to fpga/build/<top>_wl<WL>/.
# Untested: written without access to Vivado.
if {[llength $argv] < 1} {
    puts "usage: vivado -mode batch -source fpga/synth_ooc.tcl -tclargs <top> \[part\] \[clock_ns\] \[WL\]"
    exit 1
}
set top      [lindex $argv 0]
set part     [expr {[llength $argv] > 1 ? [lindex $argv 1] : "xc7a100tcsg324-1"}]
set clock_ns [expr {[llength $argv] > 2 ? [lindex $argv 2] : 10.0}]
set wl       [expr {[llength $argv] > 3 ? [lindex $argv 3] : 12}]

set root [file normalize [file dirname [info script]]/..]
set rtl  $root/rtl
set out  $root/fpga/build/${top}_wl${wl}
file mkdir $out

read_verilog [glob $rtl/*.v]

set xdc $out/clock.xdc
set fh [open $xdc w]
puts $fh "create_clock -name clk -period $clock_ns \[get_ports clk\]"
close $fh
read_xdc -mode out_of_context $xdc

synth_design -top $top -part $part -mode out_of_context -include_dirs $rtl -generic WL=$wl

write_checkpoint -force $out/${top}_synth.dcp
report_utilization    -file $out/${top}_utilization.rpt
report_timing_summary -file $out/${top}_timing.rpt
puts "reports written to $out"
