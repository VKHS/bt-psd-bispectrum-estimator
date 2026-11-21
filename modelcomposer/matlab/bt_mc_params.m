function p = bt_mc_params()
% Parameters of the generated Black Box 'bt_psd_mc' (python/scripts/gen_mc.py; do not edit, regenerate instead).
%   wl         word length of the samples: Fix_<wl>_<wl-1> at the Black Box input (-1 <= x < 1)
%   nchunks    chunks of 512 samples recorded per run;  segments = nchunks - 1 overlapping 1024-sample segments are averaged
%   avg_shift  ceil(log2(segments)): the hardware divides the accumulated power by 2^avg_shift
%   accw       width of rd_data (UFix_<accw>_0), bins = 1024 (two-sided spectrum)
%   scale      P_welch = rd_data * scale  (btpsd.estimator.hw_power_scale): the scale of the floating-point Welch average with the
%              '/L' normalisation of btpsd.reference.welch_twosided(scaling = 'paper')
p.wl = 12;
p.round = 1;
p.sq_serial = 1;
p.nchunks = 5;
p.nsamp = 2560;
p.segments = 4;
p.avg_shift = 2;
p.accw = 32;
p.bins = 1024;
p.in_type = 'Fix_12_11';
p.out_type = 'UFix_32_0';
p.scale = 0.000244140625;
end
