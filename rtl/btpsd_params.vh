// Pipeline latencies shared by the RTL (cycles from presenting the inputs to valid outputs).
// They are structural: cru_rotate has four register stages, r8_butterfly = 1 + LAT_CRU + 1 + LAT_CRU.
`ifndef BTPSD_PARAMS_VH
`define BTPSD_PARAMS_VH
`define LAT_CRU   4
`define LAT_BFLY  10
`endif
