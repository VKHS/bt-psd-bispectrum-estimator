function mdl = bt_mc_build(source, varargin)
%BT_MC_BUILD  Build the Simulink / Model Composer test bench of the PSD estimator Black Box.
%
%   mdl = bt_mc_build('sine')                    Sine Wave + Random Number sources (the default)
%   mdl = bt_mc_build('signal')                  Signal Generator ('Waveform' option: sine, square, sawtooth)
%   mdl = bt_mc_build('chirp')                   Chirp Signal, 0.01 to 0.45 cycles per sample over the record
%   mdl = bt_mc_build('file', 'File', f)         From Workspace, f = .mat written by python/scripts/export_mc_stimulus.py
%
%   Options (name, value):  'Name' model name (bt_psd_mc_demo), 'ToneBin' tone in bins of the 1024-point spectrum (100.37),
%   'ToneAmp' (0.5, fraction of full scale), 'NoiseRms' (0.05, Random Number standard deviation), 'Waveform' ('sine'),
%   'File' (stimulus file for the source 'file').
%
%   The test bench (all at one sample period = one clock cycle):
%     source -> Gateway In (Fix_WL_(WL-1), round half away from zero, saturate) -> Black Box bt_psd_mc.in_data
%     Constant 1 -> Gateway In (Boolean) -> in_valid          Constant 0 -> Gateway In (Boolean) -> rst
%     From Workspace (address ramp 0..1023 repeating) -> Gateway In (UFix_10_0) -> rd_addr
%     busy, done, rd_data -> Gateway Out -> To Workspace;  the gateway outputs of in_data, in_valid and rd_addr are logged as well
%     (bt_in, bt_valid, bt_addr, bt_rd_data, bt_done, bt_busy).  The model is saved to modelcomposer/work.
%
%   The Black Box is switched from its default simulation mode 'Inactive' (outputs zeros) to 'Vivado Simulator': the HDL is
%   co-simulated by the Vivado simulator that comes with the installation.
%
%   NOT EXECUTED where this repository was written (there is no MATLAB, Simulink or Model Composer there): the block library
%   paths and mask parameter names below are those of the AMD HDL blockset (library xbsIndex_r4) and of Simulink as documented,
%   every parameter that is set by name is set inside try/catch, and whatever could not be set or connected is listed at the
%   end under "manual steps"; docs/MODEL_COMPOSER.md section 4 gives the same bench block by block.
%
%   Start MATLAB through Model Composer first (see the guide), so that the HDL blockset is on the path.

    if nargin < 1
        source = 'sine';
    end
    ip = inputParser;
    addParameter(ip, 'Name', 'bt_psd_mc_demo');
    addParameter(ip, 'File', '');
    addParameter(ip, 'ToneBin', 100.37);
    addParameter(ip, 'ToneAmp', 0.5);
    addParameter(ip, 'NoiseRms', 0.05);
    addParameter(ip, 'Waveform', 'sine');
    parse(ip, varargin{:});
    o = ip.Results;

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    p = bt_mc_params();
    mdl = o.Name;
    notes = {};

    try
        load_system('xbsIndex_r4');
    catch err
        error('bt_mc_build:library', ['The AMD HDL blockset (library xbsIndex_r4) cannot be loaded: %s\n' ...
              'Start MATLAB through Model Composer (source Vitis/<version>/settings64.sh and Model_Composer/settings64.sh, ' ...
              'then run model_composer) and try again.'], err.message);
    end
    if bdIsLoaded(mdl)
        close_system(mdl, 0);
    end
    new_system(mdl);
    open_system(mdl);

    X = [30 200 380 560 780 1000];
    rect = @(c, r, w, h) [X(c) 40 + 70 * r X(c) + w 40 + 70 * r + h];

    % ---- data source ---------------------------------------------------------------------------------------------------
    src_out = '';
    switch lower(source)
        case 'sine'
            add_blk('simulink/Sources/Sine Wave', 'sine', rect(1, 0, 60, 40), 'SineType', 'Time based', ...
                    'Amplitude', num2str(o.ToneAmp), 'Frequency', sprintf('2*pi*%.10g/1024', o.ToneBin), ...
                    'Phase', '0', 'SampleTime', '1');
            add_blk('simulink/Sources/Random Number', 'noise', rect(1, 1, 60, 40), 'Mean', '0', ...
                    'Variance', sprintf('%.10g', o.NoiseRms^2), 'Seed', '1', 'SampleTime', '1');
            add_blk('simulink/Math Operations/Sum', 'sum', rect(2, 0, 30, 40), 'Inputs', '++');
            wire_up('sine/1', 'sum/1');
            wire_up('noise/1', 'sum/2');
            src_out = 'sum/1';
            solver = 'FixedStepDiscrete';
        case 'signal'
            add_blk('simulink/Sources/Signal Generator', 'siggen', rect(1, 0, 60, 40), 'WaveForm', o.Waveform, ...
                    'Amplitude', num2str(o.ToneAmp), 'Frequency', sprintf('%.10g/1024', o.ToneBin), 'Units', 'Hertz');
            src_out = 'siggen/1';
            solver = 'ode1';
        case 'chirp'
            add_blk('simulink/Sources/Chirp Signal', 'chirp', rect(1, 0, 60, 40), 'f1', '0.01', ...
                    'T', num2str(p.nsamp), 'f2', '0.45');
            src_out = 'chirp/1';
            solver = 'ode1';
        case 'file'
            f = o.File;
            if isempty(f)
                error('bt_mc_build:file', 'give the stimulus file: bt_mc_build(''file'', ''File'', ''modelcomposer/stimulus/<name>.mat'')');
            end
            if ~(startsWith(f, filesep) || ~isempty(regexp(f, '^[A-Za-z]:', 'once')))
                f = fullfile(pwd, f);
            end
            if ~exist(f, 'file')
                error('bt_mc_build:file', 'stimulus file %s not found; write it with python/scripts/export_mc_stimulus.py', f);
            end
            evalin('base', sprintf('load(''%s'');', f));
            set_param(mdl, 'InitFcn', sprintf('load(''%s'');', f));
            add_blk('simulink/Sources/From Workspace', 'stim', rect(1, 0, 70, 40), 'VariableName', 'bt_stim', ...
                    'SampleTime', '1', 'Interpolate', 'off', 'OutputAfterFinalValue', 'Setting to zero');
            src_out = 'stim/1';
            solver = 'FixedStepDiscrete';
        otherwise
            error('bt_mc_build:source', 'unknown source ''%s'' (sine, signal, chirp, file)', source);
    end

    % ---- gateways in and the Black Box ----------------------------------------------------------------------------------
    gw_in('gw_in_data', rect(3, 0, 70, 40), 'Signed  (2''s comp)', p.wl, p.wl - 1);
    wire_up(src_out, 'gw_in_data/1');

    add_blk('simulink/Sources/Constant', 'one', rect(1, 3, 40, 30), 'Value', '1');
    gw_in('gw_in_valid', rect(3, 3, 70, 30), 'Boolean', 1, 0);
    wire_up('one/1', 'gw_in_valid/1');

    add_blk('simulink/Sources/Constant', 'zero', rect(1, 4, 40, 30), 'Value', '0');
    gw_in('gw_in_rst', rect(3, 4, 70, 30), 'Boolean', 1, 0);
    wire_up('zero/1', 'gw_in_rst/1');

    add_blk('simulink/Sources/From Workspace', 'addr_src', rect(1, 5, 70, 40), 'VariableName', 'bt_addr_stim', ...
            'SampleTime', '1', 'Interpolate', 'off', 'OutputAfterFinalValue', 'Holding final value');
    gw_in('gw_in_addr', rect(3, 5, 70, 30), 'Unsigned', 10, 0);
    wire_up('addr_src/1', 'gw_in_addr/1');

    add_blk('xbsIndex_r4/Black Box', 'bt_psd_mc', [X(4) 40 X(4) + 170 40 + 230]);
    set_blk('bt_psd_mc', 'config_file', fullfile(here, 'bt_psd_mc_config.m'));
    % A Black Box starts in the simulation mode 'Inactive' (it ignores its inputs and outputs zeros); the HDL has to be
    % co-simulated in the Vivado simulator.  The mask parameter is looked up by its prompt, because its name is not documented.
    set_sim_mode('bt_psd_mc', 'Vivado Simulator');
    % port order of the Black Box = order of addSimulinkInport / addSimulinkOutport in bt_psd_mc_config.m
    wire_up('gw_in_rst/1', 'bt_psd_mc/1');
    wire_up('gw_in_valid/1', 'bt_psd_mc/2');
    wire_up('gw_in_data/1', 'bt_psd_mc/3');
    wire_up('gw_in_addr/1', 'bt_psd_mc/4');

    % ---- logging: gateway out -> To Workspace --------------------------------------------------------------------------
    taps = {'bt_psd_mc/1', 'busy', 'bt_busy'; 'bt_psd_mc/2', 'done', 'bt_done'; 'bt_psd_mc/3', 'rd_data', 'bt_rd_data'; ...
            'gw_in_data/1', 'in_data', 'bt_in'; 'gw_in_valid/1', 'in_valid', 'bt_valid'; 'gw_in_addr/1', 'rd_addr', 'bt_addr'};
    for k = 1:size(taps, 1)
        gname = ['gw_out_' taps{k, 2}];
        wname = ['log_' taps{k, 2}];
        add_blk('xbsIndex_r4/Gateway Out', gname, [X(5) 40 + 60 * (k - 1) X(5) + 60 40 + 60 * (k - 1) + 30]);
        add_blk('simulink/Sinks/To Workspace', wname, [X(6) 40 + 60 * (k - 1) X(6) + 90 40 + 60 * (k - 1) + 30], ...
                'VariableName', taps{k, 3}, 'SaveFormat', 'Array', 'MaxDataPoints', 'inf', 'SampleTime', '-1');
        wire_up(taps{k, 1}, [gname '/1']);
        wire_up([gname '/1'], [wname '/1']);
    end

    % ---- the hub / token of the HDL flow ------------------------------------------------------------------------------
    hub = find_system('xbsIndex_r4', 'SearchDepth', 1, 'RegExp', 'on', 'Name', 'Model Composer Hub|System Generator');
    if isempty(hub)
        notes{end + 1} = ['add the Vitis Model Composer Hub block (search for "Model Composer Hub" in the Library Browser) to the model; ' ...
                          'set its sample period (Simulink system period) to 1'];
    else
        try
            add_block(hub{1}, [mdl '/hub'], 'Position', [X(1) 40 + 70 * 6 X(1) + 80 40 + 70 * 6 + 40]);
        catch err
            notes{end + 1} = sprintf('add the Model Composer Hub block by hand (%s)', err.message);
        end
    end

    % ---- model settings ---------------------------------------------------------------------------------------------
    set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', solver, 'FixedStep', '1', 'StopTime', num2str(p.nsamp + 60000 * p.segments));
    workdir = fullfile(here, '..', 'work');
    if ~exist(workdir, 'dir')
        mkdir(workdir);
    end
    save_system(mdl, fullfile(workdir, [mdl '.slx']));

    fprintf('Built %s (source: %s), saved to %s\n', mdl, source, fullfile(workdir, [mdl '.slx']));
    fprintf('Black Box: WL %d, %d chunks (%d segments, avg_shift %d), read data %s\n', p.wl, p.nchunks, p.segments, p.avg_shift, p.out_type);
    if isempty(notes)
        fprintf('Manual steps: none were reported by the script; still check the gateways against the table in docs/MODEL_COMPOSER.md.\n');
    else
        fprintf('Manual steps (the script could not do these):\n');
        for k = 1:numel(notes)
            fprintf('  %d. %s\n', k, notes{k});
        end
    end
    fprintf('Next: bt_mc_run(''%s'')\n', mdl);

    % ======================================================================================================================
    function add_blk(src, name, position, varargin)
        try
            add_block(src, [mdl '/' name], 'Position', position);
        catch err
            error('bt_mc_build:block', 'cannot add the block %s: %s', src, err.message);
        end
        for i = 1:2:numel(varargin)
            set_blk(name, varargin{i}, varargin{i + 1});
        end
    end

    function set_blk(name, param, value)
        try
            set_param([mdl '/' name], param, value);
        catch err
            notes{end + 1} = sprintf('%s: set %s = %s by hand (%s)', name, param, value, strtrim(err.message));
        end
    end

    function gw_in(name, position, arith, nbits, binpt)
        add_blk('xbsIndex_r4/Gateway In', name, position);
        set_blk(name, 'arith_type', arith);
        if ~strcmp(arith, 'Boolean')
            set_blk(name, 'n_bits', num2str(nbits));
            set_blk(name, 'bin_pt', num2str(binpt));
        end
        if strncmp(arith, 'Signed', 6)
            set_blk(name, 'quantization', 'Round  (unbiased: +/- Inf)');
            set_blk(name, 'overflow', 'Saturate');
        end
        set_blk(name, 'period', '1');
    end

    function set_sim_mode(name, mode)
        try
            blk = [mdl '/' name];
            dp = get_param(blk, 'DialogParameters');
            fn = fieldnames(dp);
            hit = '';
            for i = 1:numel(fn)
                if ~isempty(regexpi(dp.(fn{i}).Prompt, 'simulation\s*mode', 'once'))
                    hit = fn{i};
                    break;
                end
            end
            if isempty(hit)
                error('the block has no dialog parameter with the prompt "Simulation mode"');
            end
            set_param(blk, hit, mode);
        catch err
            notes{end + 1} = sprintf(['%s: set "Simulation mode" to "%s" in the block dialog; the default "Inactive" makes the block ' ...
                                      'output zeros (%s)'], name, mode, strtrim(err.message));
        end
    end

    function wire_up(from, to)
        try
            add_line(mdl, from, to, 'autorouting', 'on');
        catch err
            notes{end + 1} = sprintf('connect %s to %s by hand (%s)', from, to, strtrim(err.message));
        end
    end
end
