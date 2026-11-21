function res = bt_mc_run(mdl, varargin)
%BT_MC_RUN  Simulate the test bench of bt_mc_build, save the logs and check the result against the Python model.
%
%   res = bt_mc_run                       simulates 'bt_psd_mc_demo'
%   res = bt_mc_run('my_model')
%
%   Options (name, value):  'StopCycles' simulation length in clock cycles (default: bt_stop_cycles of a stimulus file when it
%   is in the base workspace, else a first guess that is doubled until done has been high for the whole read-out),
%   'Check' run python/scripts/compare_mc_result.py afterwards (true), 'Python' interpreter (python3).
%
%   The address input of the Black Box is a ramp 0..1023 that repeats every 1024 clock cycles (From Workspace bt_addr_stim, written
%   here for the chosen length).  While done = 1 the ramp reads every bin; rd_data follows rd_addr one clock cycle later.  The logs
%   (one entry per clock cycle) are saved as the struct bt_log in modelcomposer/work/bt_log.mat:
%       in_code    integer code at the Gateway In output (value * 2^(WL-1))      in_valid   valid input of the block
%       rd_addr    address input        rd_data   result word        done, busy   outputs of the block        wl, nchunks
%   The expected spectrum is computed by the checker from in_code, so it does not depend on the source block.
%
%   NOT EXECUTED where this repository was written (no MATLAB there); see docs/MODEL_COMPOSER.md.

    if nargin < 1
        mdl = 'bt_psd_mc_demo';
    end
    ip = inputParser;
    addParameter(ip, 'StopCycles', 0);
    addParameter(ip, 'Check', true);
    addParameter(ip, 'Python', 'python3');
    parse(ip, varargin{:});

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    root = fullfile(here, '..', '..');
    workdir = fullfile(here, '..', 'work');
    p = bt_mc_params();

    if ~bdIsLoaded(mdl)
        slx = fullfile(workdir, [mdl '.slx']);
        if exist(slx, 'file')
            load_system(slx);
        else
            load_system(mdl);
        end
    end

    stop = ip.Results.StopCycles;
    if stop <= 0
        if evalin('base', 'exist(''bt_stop_cycles'', ''var'')')
            stop = evalin('base', 'bt_stop_cycles');
        else
            stop = p.nsamp + p.segments * 45000 + 4000;
        end
    end

    ok = false;
    for attempt = 1:6
        t = (0:stop);
        t = t(:);
        assignin('base', 'bt_addr_stim', [t mod(t, 1024)]);
        set_param(mdl, 'StopTime', num2str(stop));
        fprintf('Simulating %s for %d clock cycles (attempt %d) ...\n', mdl, stop, attempt);
        sim(mdl);
        done = evalin('base', 'bt_done');
        if nnz(done ~= 0) >= 1100
            ok = true;
            break;
        end
        if any(done ~= 0)
            fprintf('done rose only %d cycles before the end; the read-out needs 1025. Doubling the length.\n', nnz(done ~= 0));
        else
            fprintf('done did not rise within %d cycles. Doubling the length.\n', stop);
        end
        stop = 2 * stop;
    end
    if ~ok
        error('bt_mc_run:done', ['done never stayed high long enough. Check that the Black Box simulates (simulation mode of the ' ...
              'block), that rst is 0 and in_valid is 1, and that the source is not constant zero.']);
    end

    from_base = @(name) evalin('base', name);
    res = struct();
    res.in_code = round(from_base('bt_in') * 2^(p.wl - 1));
    res.in_valid = from_base('bt_valid');
    res.rd_addr = from_base('bt_addr');
    res.rd_data = from_base('bt_rd_data');
    res.done = from_base('bt_done');
    res.busy = from_base('bt_busy');
    res.wl = p.wl;
    res.nchunks = p.nchunks;

    if ~exist(workdir, 'dir')
        mkdir(workdir);
    end
    matfile = fullfile(workdir, 'bt_log.mat');
    bt_log = res;
    save(matfile, 'bt_log', '-v7');
    fprintf('Saved %s (%d clock cycles)\n', matfile, numel(res.done));

    if ip.Results.Check
        script = fullfile(root, 'python', 'scripts', 'compare_mc_result.py');
        setenv('PYTHONPATH', fullfile(root, 'python'));
        cmd = sprintf('%s "%s" "%s"', ip.Results.Python, script, matfile);
        if isunix
            cmd = ['LD_LIBRARY_PATH= ' cmd];
        end
        [status, out] = system(cmd);
        fprintf('%s\n', out);
        if status ~= 0
            fprintf(['The checker returned %d (1: the result differs from the model or the log is unusable, 2: the file could ' ...
                     'not be read, 127: no Python found). Run it by hand from a shell:\n  PYTHONPATH=python python3 %s %s\n'], ...
                    status, 'python/scripts/compare_mc_result.py', 'modelcomposer/work/bt_log.mat');
        end
    end
end
