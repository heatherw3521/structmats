function mp_memo1_run(what, opts)
%MP_MEMO1_RUN  Run the experiments of memo 1 and write their results as text.
%   Run from hss_examples/minnorm_paper:
%       mp_memo1_run('timing')     timing and large-n accuracy (Section 5, Figures 1, 3)
%                                  -> results/memo1_scaling.txt
%       mp_memo1_run('accuracy')   Table 2 and Figure 2
%                                  -> results/memo1_correctness.txt, memo1_conditioning.txt
%       mp_memo1_run               both
%   mp_memo1_run(what, opts) passes opts to the experiments, e.g.
%   opts.maxexp = 16 (largest n = 2^16 in the timing run; default 17),
%   opts.quick = true (small smoke test, *_quick.txt files).
%   The text files are plain CSV with two '#' header lines (MATLAB version,
%   platform, thread count, date). For timings: close other programs first.
%   Expected run time in MATLAB: accuracy a few minutes; timing roughly
%   10-30 minutes at the default sizes.
if nargin < 1 || isempty(what), what = 'all'; end
if nargin < 2, opts = struct(); end
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
addpath(fullfile(here, 'lib'));
if ~exist('OCTAVE_VERSION', 'builtin')
    addpath(fileparts(fileparts(here)));     % repo root: @hss, +hssutil
end
if any(strcmp(what, {'all', 'accuracy'}))
    fprintf('==== Table 2: correctness\n');  mp_exp_correctness(opts);
    fprintf('==== Figure 2: conditioning\n'); mp_exp_conditioning(opts);
end
if any(strcmp(what, {'all', 'timing'}))
    fprintf('==== Figures 1 and 3: timing and large-n accuracy\n'); mp_exp_scaling(opts);
end
fprintf('done: results in %s\n', fullfile(here, 'results'));
end
