# python_reference: independent implementation used for the memos

A Python port of `@hss/private/hss_ulvminnormsolve.m` plus the test-matrix families, the weighted
and Tikhonov extensions, and every experiment behind the figures of
`HSS_Min_Norm_Solve_Memo_v2.pdf` (memo 1) and `HSS_Weighted_Tikhonov_Memo.pdf` (memo 2).
It was cross-checked against the repo's MATLAB solver, current and legacy, on identical matrices
(memo 1, Table 6).

Requirements: Python 3.11, NumPy, SciPy, matplotlib, mpmath (60-digit references only).

## Package `hssmn/`

| file | contents |
|---|---|
| `hss.py` | generator-form HSS matrix (`HSS`), products, adjoint, densify, leaf merge, constructor tree rule (`matlab_tree`) |
| `ulv.py` | the solver: `MinNormFactor(H, slack='relax'|'collapse', root='factored'|'pinv')`, `.solve(b)` |
| `build.py` | deterministic IDs, constructor rule (`mode='full'`), O(n) proxy builder (`mode='proxy'`), F1 `lrbd_hss`, F2 `random_hss` |
| `families.py` | F3 Cauchy, F4 convolution, F5 NUDFT Cauchy-like, F6 Toeplitz to Cauchy-like |
| `variants.py` | memo 2: `weighted_minnorm`, `augment`, `tikhonov` |
| `refs.py` | dense QR / SVD min-norm, normal equations, CGNE, 60-digit `mpmath` reference, error measures |
| `cost.py` | flop count of one factorization (memo 1, Sec. 5) |
| `export.py` | export generators to `.mat` for the MATLAB cross-check |

The defaults (relaxed width, factored backward-stable root solve, economic size-reduction QR) match
the current MATLAB solver. Two options reproduce the legacy file `@hss/legacy/hss_ulvminnormsolve.m`:
`slack='collapse'` its level-collapse rule, `root='pinv'` its root solve `pinv(H.D)*b` (explicit
pseudo-inverse).

## Experiments `exp/`

Run from `exp/` with single-threaded BLAS, on an otherwise idle machine, e.g.
`OPENBLAS_NUM_THREADS=1 python3 exp_m3_scaling.py`. Results go to `exp/results/*.json`.

| script | memo | output |
|---|---|---|
| `exp_m1_correctness.py` | M1 Table 4 | `m1_correctness.json` |
| `exp_m2_conditioning.py`, `exp_m2b_deep.py` | M1 Fig. 1 | `m2_conditioning.json`, `m2b_deep.json` |
| `exp_m3_scaling.py` | M1 Figs. 2-3, Table 5 | `m3_scaling.json` |
| `exp_m5_sweeps.py`, `exp_m5a_rank.py` | M1 Fig. 4 | `m5_sweeps.json` (the second replaces the rank sweep) |
| `exp_m6_square.py` | M1 Sec. 6.7 (square case) | printed |
| `exp_a1_nudft.py`, `exp_a2_a3_a4.py` | M1 Sec. 7 | `a1_nudft.json`, `a2_a3_a4_*.json` |
| `exp_w_t.py w1 / w2 / t1 / t2 t3` | M2 Sec. 5 | `wt_*.json` |
| `exp_apps2.py bp / db / mri` | M2 Sec. 6 | `apps2_*.json` |

Figures: `python3 fig_memo1.py`, `python3 fig_memo2.py` (to `exp/figures/`). Tables:
`python3 tables.py` (LaTeX fragments for memo 1; the cross-check table reads the Octave outputs
`results/xval2_results_current.txt` (current MATLAB solver) and `results/xval2_results_legacy.txt`
(the legacy file), paths set at the top of `tables.py`).
