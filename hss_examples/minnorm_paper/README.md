# minnorm_paper — test suite for the HSS minimum-norm solver

Paper-style tests for `@hss/private/hss_ulvminnormsolve.m` (called through `H \ b` for square and
wide `H`), and for the weighted minimum-norm and Tikhonov methods of memo 2 (`minnorm`, `tikhonov`).

Companion documents (repo root):
- `HSS_Min_Norm_Solve_Memo_v2.pdf`: algorithm, proofs, complexity, numerical study (memo 1)
- `HSS_Weighted_Tikhonov_Memo.pdf`: weighted minimum norm and Tikhonov (memo 2)

## Memo 1: results as text files

```matlab
cd hss_examples/minnorm_paper
mp_memo1_run('timing')     % memo 1, Section 5, Figures 1 and 3 -> results/memo1_scaling.txt
mp_memo1_run('accuracy')   % memo 1, Table 2 and Figure 2       -> results/memo1_correctness.txt,
                           %                                       results/memo1_conditioning.txt
mp_memo1_run               % both
```
The files are plain CSV with two `#` header lines (MATLAB version, platform, threads, date).
`notes/minnorm_v2/make_figures.py` turns them into the memo's figures, Table 2 and the numbers
quoted in the text. Timing: close other programs first; the timing run takes roughly 10-30 minutes
at the default sizes (n up to 2^17; `mp_memo1_run('timing', struct('maxexp', 16))` for less).

## Run

```matlab
cd hss_examples/minnorm_paper
mp_selftest                       % 33 unit checks, small sizes (~1 min)
mp_test_solver                    % 49 regression checks: H\b, stored factors, minnorm, tikhonov (~2 min)
opts.quick = true;  run_minnorm_suite   % smoke test of every experiment (minutes)
clear opts; opts.maxexp = 17; opts.reps = 3; run_minnorm_suite   % full run
```
Results go to `results/*.mat`, figures to `figures/*.pdf|png` (`mp_plot_results` re-plots
from saved results). Each experiment is also callable alone, e.g. `S = mp_exp_scaling(opts)`.
Timings of `H\b` are taken with `lib/mp_time_solve.m`, which separates the first solve (factor and
solve) from a repeat solve that uses the factors `H` keeps.

## What each experiment does

| file | memo section | what it measures |
|---|---|---|
| `mp_test_solver.m` | -- | regression checks of `H\b`: every family vs dense min-norm (error, backward error, no merged levels), square, complex, several right-hand sides, stored factors (reuse, `clearfactors`, dropped on modification), `minnorm`/`tikhonov` with vector and leaf-block weights vs dense formulas, named errors |
| `mp_exp_correctness.m` | M1 7.3, Table 2 | every family vs dense min-norm (residual, solver error, null-space fraction, error vs the original matrix, κ, slack condition) |
| `mp_exp_conditioning.m` | M1 7.4, Fig. 2 | accuracy vs κ (graded F2, Gaussian blur): ULV vs dense QR, normal equations, CGNE |
| `mp_exp_scaling.m` | M1 5 and 7.5, Figs. 1 and 3 | first and repeat solve times and manufactured-solution errors to n = 2^maxexp (default leaf size 200); dense QR and CGNE on F2 |
| `mp_exp_sweeps.m` | -- | time vs rank, blocksize, aspect ratio (not in memo 1) |
| `mp_exp_applications.m` | -- | NUDFT interpolation, deconvolution, Toeplitz→Cauchy-like, 2-D blur (not in memo 1) |
| `mp_exp_weighted_tikhonov.m` | M2 5 | W1/W2 weighted accuracy and scale, T1 λ-path accuracy (wide/square/tall, general form), T2 scale |
| `mp_exp_apps2.m` | M2 6 | basis pursuit (IRLS, Douglas–Rachford), Tikhonov deblurring, k-space-weighted MRI |

## Test matrices (exact definitions in memo 1, Table 1)

| family | generator | dense ever formed? |
|---|---|---|
| F1 LR + block diagonal (showcase `rectsampler`) | `lib/mp_gen_lrbd` | no (generator form, O(nk)) |
| F2 random HSS, rank k, real/complex | `lib/mp_gen_random` | no |
| F3 interlaced Cauchy 1/(x_i − y_j) | `lib/mp_kernel_cauchy` + `mp_hss_kernel` | no (O(n) proxy build) |
| F4 decimated Gaussian/Ricker convolution | `lib/mp_kernel_conv` + `mp_hss_kernel` | no (banded) |
| F5 NUDFT Cauchy-like (Dirichlet kernel) | `lib/mp_kernel_nudft` + `mp_hss_kernel` | no |
| F6 Toeplitz → Cauchy-like (unitary FFT transform) | `lib/mp_kernel_toeplitz` + `mp_hss_kernel` | no |

Dense matrices appear only as *references* at moderate size (n ≤ 2^13). Large-n accuracy uses
manufactured solutions: `z ~ N(0,I)`, `x* = H'*z`, `b = H*x*`. Then `x*` is exactly the
minimum-norm solution of the HSS system (`lib/mp_manufactured`).

## Library (`lib/`)

- `mp_hss_from_generators` / `mp_hss_to_generators`: build an `@hss` object from generators
  (D, U/R, V/W, B) and back. The layout is identical to `hss_constructor.m` output, so `H*x`,
  `H'`, `H\b` and `full(H)` work unchanged.
- `mp_tree`: the constructor's tree rule. `mp_tree_aligned`: geometry-aligned tree for
  nonuniform sampling (used by `mp_exp_applications`).
- `mp_hss_kernel`: nested-ID construction from entry evaluation. `mode='proxy'` (default) is
  O(n); `mode='full'` is the constructor's O(mn) rule.
- `mp_id_rows`, `mp_id_cols`, `mp_pivqr`: deterministic IDs (same relative-threshold rule as
  `inter_decompv3`, no sketch, no rank cap).
- Memo 2: `mp_weighted_minnorm`, `mp_tikhonov` (thin wrappers of the `@hss` methods `minnorm` and
  `tikhonov`), and the generator edits `mp_gen_scale`, `mp_gen_rightblock`, `mp_gen_augment` (used by
  the experiments; the methods have their own copies in `@hss/private`).
- Helpers: `mp_minnorm_dense` (QR reference), `mp_manufactured`, `mp_slack_ok`, `mp_timeit`,
  `mp_rel`, `mp_nullfrac`, `mp_save`, `mp_figsave`, `mp_style`, `mp_conv_apply`.

## Solver versions

- `@hss/private/hss_ulvminnormsolve.m` (used by `H\b` for square and wide `H`): the current solver.
  It has five changes over the legacy version: backward-stable root solve (SVD factors
  instead of `pinv(H.D)*b`), rank checks with named errors, relaxed width `p' = min(p, l+n)` (no slack
  precondition, no level merging), economic QR in the size reduction, and a factor/solve split.
  Leaves where the size reduction would discard nothing skip it, which makes it the square solver too.
  Several right-hand sides can be passed as columns of `b`.
- Stored factors: `H` keeps its factorization after the first `H\b`; every later `H\b` with the
  same `H` only solves, O(n r) (repeated projections in Douglas-Rachford/ADMM, new right-hand sides).
  Nothing to call. Assigning to any part of `H` drops them (a copy `H2 = H` shares them until one is
  modified); they are not saved to `.mat` files; `clearfactors(H)` frees the memory.
- `x = minnorm(H, b, 'Weight', L)` and `[x, r] = tikhonov(H, b, lambda, 'Weight', L, 'DataWeight', S)`
  (memo 2): `L`, `S` are vectors or cell arrays of leaf blocks. They keep the factors for the last
  weight / `(lambda, L, S)`, so repeat calls with the same values only solve.
- `@hss/legacy/`: `hss_ulvminnormsolve.m`, the earlier version (explicit `pinv` at the root,
  slack check + level merging); `hss_ulvvecsolve.m`, the earlier square solver (one right-hand side,
  no stored factors); `minnormfactor.m`, superseded by the stored factors. To run the first, copy it to a folder on the MATLAB path outside
  `@hss` and call `hss_ulvminnormsolve(H, b)` from there (MATLAB does not search subfolders of a
  class folder); `H\b` always uses the private (current) version.
- `proposed/mp_ulvminnormsolve_relaxed.m`: the standalone variant written before the update
  (relaxed width and factored root only, full QR, no factor/solve split). Kept for reference; it
  agrees with `H\b` to rounding.

## Validation status

No MATLAB licence was available while this was written. Every file here was executed under
GNU Octave 8.4 against the repo's own `@hss` code (with a `dictionary` shim and an
Octave-parsable copy of the `arguments` blocks; see `octave_compat/`): `mp_selftest` passes
33/33, `mp_test_solver` passes 49/49, and the quick suite (`opts.quick = true`) runs end to end, including every figure.
Please run `mp_selftest` once in MATLAB first.

Notes on reading the output:
- CGNE uses `pcg`; when `flag ~= 0` it did not converge, and `pcg` returns its minimal-residual
  iterate (often the starting guess for very ill-conditioned problems, hence `fwd 1.0`).
- Octave timings are not representative (its HSS products are much slower than MATLAB's); the
  timings in memo 1 come from MATLAB (`mp_memo1_run('timing')`, results/memo1_scaling.txt).
