"""E2: stability vs conditioning.  Forward error against a 60-digit reference
and normwise backward error, for the HSS ULV solve, dense QR / SVD min-norm,
dense normal equations (Cholesky of HH^*) and CG on HH^* (CGNE, HSS matvecs).

Families (exact HSS, so all solvers see the same matrix):
  (a) F2 random HSS times a random column grading 10^{-alpha u_j}, u_j ~ U[0,1]
  (b) F4 'valid' Gaussian convolution (stride 1) with growing blur width sigma
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from hssmn import *
from hssmn.families import *
from hssmn.refs import *
from hssmn.hss import matlab_tree
from common import dump

rng = np.random.default_rng(7)
U_ROUND = np.finfo(float).eps / 2


def solve_all(H, b, xref):
    A = H.dense()
    nA = np.linalg.norm(A, 2)
    res = {}
    def rec(key, x, extra=None):
        res[key] = dict(fwd=rel(x, xref), bwd=float(np.linalg.norm(b - A @ x) / (nA * np.linalg.norm(x) + np.linalg.norm(b))))
        if extra:
            res[key].update(extra)
    rec('ulv', MinNormFactor(H, slack='relax').solve(b))
    rec('ulv_pinv', MinNormFactor(H, slack='relax', root='pinv').solve(b))   # repo root: pinv(D)*b
    rec('qr', minnorm_qr(A, b))
    rec('svd', minnorm_pinv(A, b))
    try:
        rec('normal_eq', minnorm_normal_eq(A, b))
    except np.linalg.LinAlgError:
        res['normal_eq'] = dict(fwd=np.nan, bwd=np.nan, failed=True)
    x, it, info = minnorm_cgne(H, b, rtol=1e-14, maxiter=2000)
    rec('cgne', x, dict(iters=it, info=int(info)))
    return res


out = dict(graded=[], blur=[])
m, n, k, bs = 128, 256, 6, 16
H0 = random_hss(m, n, k, bs, rng)
u = rng.random(n)
for alpha in np.arange(0, 26, 2.0):
    H = H0.scale_columns(10.0 ** (-alpha * u))
    A = H.dense()
    s = np.linalg.svd(A, compute_uv=False)
    b = rng.standard_normal(m)
    xref = minnorm_mp(A, b, dps=60)
    r = solve_all(H, b, xref)
    r.update(alpha=float(alpha), kappa=float(s[0] / s[-1]))
    out['graded'].append(r)
    print('graded alpha=%4.1f kappa=%.2e  ulv %.1e (bwd %.1e, pinv-root bwd %.1e)  qr %.1e  ne %.1e  cgne %.1e (%d its)' %
          (alpha, r['kappa'], r['ulv']['fwd'], r['ulv']['bwd'], r['ulv_pinv']['bwd'], r['qr']['fwd'], r['normal_eq']['fwd'], r['cgne']['fwd'], r['cgne']['iters']), flush=True)

for sig in [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4]:
    sp = conv_spec(300, stride=1, kind='gauss', sigma=sig, halfwidth=int(np.ceil(5 * sig)))
    L, rb, cb = matlab_tree(sp.m, sp.n, 16)
    H = hss_build(sp, L, rb, cb, 1e-15, mode='proxy')
    A = H.dense()
    s = np.linalg.svd(A, compute_uv=False)
    b = rng.standard_normal(sp.m)
    xref = minnorm_mp(A, b, dps=60)
    r = solve_all(H, b, xref)
    r.update(sigma=sig, kappa=float(s[0] / s[-1]), m=sp.m, n=sp.n)
    out['blur'].append(r)
    print('blur sigma=%.1f %dx%d kappa=%.2e  ulv %.1e (bwd %.1e, pinv-root bwd %.1e)  qr %.1e  ne %.1e  cgne %.1e (%d its)' %
          (sig, sp.m, sp.n, r['kappa'], r['ulv']['fwd'], r['ulv']['bwd'], r['ulv_pinv']['bwd'], r['qr']['fwd'], r['normal_eq']['fwd'], r['cgne']['fwd'], r['cgne']['iters']), flush=True)
out['u'] = U_ROUND
dump('m2_conditioning', out)
