"""A1: band-limited (min-energy) interpolation from nonuniform samples.

Model: f(x) = sum_{k=-N/2}^{N/2-1} c_k e^{2 pi i k x}, samples f(x_j), j=1..m < N.
A_jk = e^{2 pi i k x_j};  C = A F^*  (F unitary DFT) is Cauchy-like and HSS;
min ||c|| s.t. Ac = f  <=>  min ||y|| s.t. C y = f,  y = F c = sqrt(1/N)*(f on the grid l/N).
Seismic reading: x_j = irregular trace positions, y = regularised traces.

Also: index-based (MATLAB constructor) tree vs geometry-aligned tree when the
sampling has a gap (missing traces).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from hssmn import *
from hssmn.families import *
from hssmn.refs import rel, minnorm_qr
from hssmn.hss import matlab_tree, tree_from_leaves
from common import dump

rng = np.random.default_rng(11)
N = 1024
k = np.arange(-N // 2, N // 2)
# real band-limited test signal with a decaying spectrum
c = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / (1 + np.abs(k) / 40.0) ** 1.5
# index i <-> k = i - N/2 ; enforce c_{-k} = conj(c_k) (real signal), drop unpaired k = -N/2
c[N // 2 + 1:] = np.conj(c[1:N // 2][::-1])
c[N // 2] = c[N // 2].real
c[0] = 0.0
f = lambda x: np.real(np.exp(2j * np.pi * np.outer(x, k)) @ c)

out = {}
for case in ['jitter', 'gap']:
    m = int(0.7 * N)
    xs = nudft_points(m if case == 'jitter' else int(m * 1.15), rng)
    if case == 'gap':
        xs = xs[(xs < 0.42) | (xs > 0.52)][:m]     # 10% of the interval has no samples
        m = len(xs)
    sp = nudft_spec(xs, N)
    fx = f(xs).astype(complex)
    # (i) MATLAB-style index tree;  (ii) geometry-aligned tree (rows split where the columns split)
    L, rb, cb = matlab_tree(m, N, 32)
    H_idx = hss_build(sp, L, rb, cb, 1e-12, mode='proxy')
    grid = np.arange(N) / N
    rbL = np.searchsorted(xs, grid[cb[L][:-1]]); rbL = np.append(rbL, m); rbL[0] = 0
    Lg, rbg, cbg = tree_from_leaves(rbL, cb[L])
    H_geo = hss_build(sp, Lg, rbg, cbg, 1e-12, mode='proxy')
    C = sp.dense()
    y_ref = minnorm_qr(C, fx)
    res = {}
    for nm, H in [('index_tree', H_idx), ('geometry_tree', H_geo)]:
        F = MinNormFactor(H, slack='relax')
        y = F.solve(fx)
        res[nm] = dict(maxrank=H.max_rank(), storage=H.storage(), err=rel(y, y_ref), res=rel(C @ y, fx),
                       collapses_matlab=int(sum(MinNormFactor(H, slack='collapse').collapses)))
        print(case, nm, res[nm], flush=True)
    y = MinNormFactor(H_idx).solve(fx)
    cmn = unitary_dft(N).conj().T @ y          # Fourier coefficients of the min-norm interpolant
    # fine-grid reconstruction for the figure
    tt = np.linspace(0, 1, 4000, endpoint=False)
    fmn = np.real(np.exp(2j * np.pi * np.outer(tt, k)) @ cmn)
    out[case] = dict(m=m, N=N, xs=xs, fx=fx.real, tt=tt, ftrue=f(tt), fmn=fmn,
                     norm_c_true=float(np.linalg.norm(c)), norm_c_mn=float(np.linalg.norm(cmn)),
                     interp_err=float(np.max(np.abs(np.real(np.exp(2j * np.pi * np.outer(xs, k)) @ cmn) - fx.real))),
                     trees=res)
    print(case, 'm=%d ||c_true||=%.3f ||c_mn||=%.3f  max|f_mn(x_j)-f(x_j)|=%.1e' %
          (m, out[case]['norm_c_true'], out[case]['norm_c_mn'], out[case]['interp_err']), flush=True)
dump('a1_nudft', out)
