# minnorm_paper — experiments for the HSS minimum-norm solver

Accuracy, conditioning, scaling and application experiments for `H \ b` (square and wide `H`,
minimum-norm solution for wide `H`) and for the weighted minimum-norm and Tikhonov methods
`minnorm` and `tikhonov`.

## Run

```matlab
cd hss_examples/minnorm_paper
mp_selftest                       % unit checks of the helpers, small sizes (~1 min)
opts.quick = true;  run_minnorm_suite   % smoke test of every experiment (minutes)
clear opts; opts.maxexp = 17; opts.reps = 3; run_minnorm_suite   % full run
```
Solver regression tests live in `hss_examples/tests/test_hss.m`
(`runtests('test_hss', 'Tag', 'families')`).

Results go to `results/*.mat`, figures to `figures/*.pdf|png` (`mp_plot_results` re-plots
from saved results). Each experiment is also callable alone, e.g. `S = mp_exp_scaling(opts)`.
Timings of `H\b` are taken with `mp_time_solve` (in `hss_examples/lib`), which separates the first
solve (factor and solve) from a repeat solve that uses the factors `H` keeps.

### Results as text

```matlab
mp_run_text_results('accuracy')   % E1, E2 -> results/correctness.txt, results/conditioning.txt
mp_run_text_results('timing')     % E3/E4  -> results/scaling.txt
mp_run_text_results               % both
```
The files are plain CSV with two `#` header lines (MATLAB version, platform, threads, date).
For timings, close other programs first; the timing run takes roughly 10-30 minutes at the default
sizes (n up to 2^17; `mp_run_text_results('timing', struct('maxexp', 16))` for less).

## Experiments

| file | what it measures |
|---|---|
| `mp_exp_correctness.m` (E1) | every family vs dense min-norm (residual, solver error, null-space fraction, error vs the original matrix, κ, leaf slack condition) |
| `mp_exp_conditioning.m` (E2) | accuracy vs κ (graded F2, Gaussian blur): ULV vs dense QR, normal equations, CGNE |
| `mp_exp_scaling.m` (E3/E4) | first and repeat solve times and manufactured-solution errors to n = 2^maxexp (default leaf size 200); dense QR and CGNE on F2 |
| `mp_exp_sweeps.m` (E5/E6) | time vs rank, blocksize, aspect ratio |
| `mp_exp_applications.m` (A1-A4) | NUDFT interpolation, deconvolution, Toeplitz→Cauchy-like, 2-D blur |
| `mp_exp_weighted_tikhonov.m` (W, T) | weighted accuracy and scale, Tikhonov λ-path accuracy (wide/square/tall, general form) and scale |
| `mp_exp_apps2.m` | basis pursuit (IRLS, Douglas–Rachford), Tikhonov deblurring, k-space-weighted MRI |

## Test matrices

| family | generator | dense ever formed? |
|---|---|---|
| F1 low rank + block diagonal | `mp_gen_lrbd` | no (generator form, O(nk)) |
| F2 random HSS, rank k, real/complex | `mp_gen_random` | no |
| F3 interlaced Cauchy 1/(x_i − y_j) | `mp_kernel_cauchy` + `mp_hss_kernel` | no (O(n) proxy build) |
| F4 decimated Gaussian/Ricker convolution | `mp_kernel_conv` + `mp_hss_kernel` | no (banded) |
| F5 NUDFT Cauchy-like (Dirichlet kernel) | `mp_kernel_nudft` + `mp_hss_kernel` | no |
| F6 Toeplitz → Cauchy-like (unitary FFT transform) | `mp_kernel_toeplitz` + `mp_hss_kernel` | no |

Dense matrices appear only as *references* at moderate size (n ≤ 2^13). Large-n accuracy uses
manufactured solutions: `z ~ N(0,I)`, `x* = H'*z`, `b = H*x*`. Then `x*` is exactly the
minimum-norm solution of the HSS system (`mp_manufactured`).

## Helpers

Generators, kernels, trees, IDs, dense references and timing are in the shared library
`hss_examples/lib` (see `help lib`). This folder's `lib/` holds only the output helpers
`mp_save`, `mp_figsave`, `mp_style`, `mp_write_table`, and `mp_slack_ok`. The scripts add both
folders to the path.

## Octave

The suite also runs under GNU Octave 8.4 or later; see `octave_compat/README.md`. Octave timings
are not representative of MATLAB's.

Reading the output: CGNE uses `pcg`; when `flag ~= 0` it did not converge, and `pcg` returns its
minimal-residual iterate (often the starting guess for very ill-conditioned problems, hence
`fwd 1.0`).
