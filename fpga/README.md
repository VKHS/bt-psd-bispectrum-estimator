# FPGA scripts

`synth_ooc.tcl` runs Vivado out-of-context synthesis of one top module and writes the utilisation and timing
summary reports to `fpga/build/` (ignored by git).

```bash
vivado -mode batch -source fpga/synth_ooc.tcl -tclargs bt_psd_top xc7a100tcsg324-1 10.0 12
vivado -mode batch -source fpga/synth_ooc.tcl -tclargs bt_bispectrum_top
vivado -mode batch -source fpga/synth_ooc.tcl -tclargs fft512_mem xc7a100tcsg324-1 5.0 16
```

* The script is untested (written without access to Vivado). The RTL passes its six testbenches in Icarus Verilog
  (`make -C sim all`); it has not been run through Vivado's simulator.
* No synthesis result exists for this repository. The resource, frequency and power numbers of the paper (Table 9) belong
  to an overlapped, multiplier-free design; the reference RTL here is sequential and uses `*` operators, so its numbers
  are not comparable (see `docs/ARCHITECTURE.md`, section 6).
* The part is a parameter; the paper names only "Artix-7".
