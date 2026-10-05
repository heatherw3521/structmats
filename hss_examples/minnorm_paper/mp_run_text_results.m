function mp_run_text_results(what, opts)
%MP_RUN_TEXT_RESULTS  Run the main experiments and write their results as text.
%   Run from hss_examples/minnorm_paper:
%       mp_run_text_results('timing')     timing and large-n accuracy (E3/E4)
%                                         -> results/scaling.txt
%       mp_run_text_results('accuracy')   correctness and conditioning (E1, E2)
%                                         -> results/correctness.txt, conditioning.txt
%       mp_run_text_results               both
%   mp_run_text_results(what, opts) passes opts to the experiments, e.g.
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
addpath(fullfile(fileparts(here), 'lib'));      % shared HSS test library (hss_examples/lib)
if ~exist('OCTAVE_VERSION', 'builtin')
    addpath(fileparts(fileparts(here)));     % repo root: @hss, +hssutil
end
if any(strcmp(what, {'all', 'accuracy'}))
    fprintf('==== E1 correctness\n');  mp_exp_correctness(opts);
    fprintf('==== E2 conditioning\n'); mp_exp_conditioning(opts);
end
if any(strcmp(what, {'all', 'timing'}))
    fprintf('==== E3/E4 timing and large-n accuracy\n'); mp_exp_scaling(opts);
end
fprintf('done: results in %s\n', fullfile(here, 'results'));
end
