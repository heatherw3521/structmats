%RUN_MINNORM_SUITE  Test suite for the HSS minimum-norm solver (H\b) and
%   the weighted minimum-norm and Tikhonov methods (minnorm, tikhonov).
%
%   Usage (from this folder):
%       opts.quick  = false;   % true: small smoke test (a few minutes)
%       opts.maxexp = 17;      % largest n = 2^maxexp in the scaling studies
%       opts.reps   = 3;       % timing repetitions (median after one warm-up)
%       run_minnorm_suite
%   Results -> results/*.mat, figures -> figures/*.pdf|png (mp_plot_results).
%
%   Experiments:
%     mp_selftest                 unit tests of the helpers
%     mp_exp_correctness   (E1)   all families vs dense min-norm
%     mp_exp_conditioning  (E2)   accuracy vs kappa
%     mp_exp_scaling       (E3/4) timing + manufactured solutions to 2^maxexp
%     mp_exp_sweeps        (E5/6) rank, blocksize, aspect ratio
%     mp_exp_applications  (A1-4) NUDFT, deconvolution, Toeplitz, 2-D
%     mp_exp_weighted_tikhonov    W1, W2, T1, T2: weighted min-norm, Tikhonov
%     mp_exp_apps2                IRLS / Douglas-Rachford, deblurring, MRI
%   mp_run_text_results runs E1-E4 alone and writes CSV text results.
%   Solver regression tests: runtests('test_hss', 'Tag', 'families').
%   No experiment forms a dense matrix above n = 2^13 (dense references are
%   only used at moderate sizes); large cases use generator-form HSS matrices
%   (mp_gen_*), or the O(n) proxy-point builder hss_from_kernel.
if ~exist('opts', 'var'), opts = struct(); end
if ~isfield(opts, 'quick'), opts.quick = false; end
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
addpath(fullfile(here, 'lib'));
addpath(fullfile(fileparts(here), 'lib'));      % shared HSS test library (hss_examples/lib)
if ~exist('OCTAVE_VERSION', 'builtin')
    addpath(fileparts(fileparts(here)));     % repo root: @hss, +hssutil
end
resdir = fullfile(here, 'results'); figdir = fullfile(here, 'figures');
tag = ''; if opts.quick, tag = '_quick'; end
fprintf('==== self test\n'); mp_selftest;
fprintf('==== E1 correctness\n');      R = mp_exp_correctness(opts);        mp_save(resdir, ['E1_correctness' tag], struct('R', R));
fprintf('==== E2 conditioning\n');     S = mp_exp_conditioning(opts);       mp_save(resdir, ['E2_conditioning' tag], S);
fprintf('==== E3/E4 scaling\n');       S = mp_exp_scaling(opts);            mp_save(resdir, ['E3_scaling' tag], S);
fprintf('==== E5/E6 sweeps\n');        S = mp_exp_sweeps(opts);             mp_save(resdir, ['E5_sweeps' tag], S);
fprintf('==== A1-A4 applications\n');  S = mp_exp_applications(opts);       mp_save(resdir, ['A_applications' tag], S);
fprintf('==== W/T weighted, Tikhonov\n'); S = mp_exp_weighted_tikhonov(opts);  mp_save(resdir, ['W_T' tag], S);
fprintf('==== apps2 applications\n');   S = mp_exp_apps2(opts);            mp_save(resdir, ['apps2' tag], S);
mp_plot_results(resdir, figdir, tag);
fprintf('done: results in %s, figures in %s\n', resdir, figdir);
