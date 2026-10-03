"""E2b: stability on a deeper tree.  F2 random HSS (2048 x 4096, rank 8,
blocksize 32, depth 6) times the column grading 10^{-alpha u_j}.  Reference:
dense QR minimum-norm solution (backward stable; a 60-digit reference is out
of reach at this size), so forward errors are meaningful down to ~kappa*u.
Both root solves are run: factored SVD (default) and the repo's explicit pinv."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from hssmn import *
from hssmn.refs import *
from common import dump

rng = np.random.default_rng(11)
m, n, k, bs = 2048, 4096, 8, 32
H0 = random_hss(m, n, k, bs, rng)
u = rng.random(n)
out = dict(rows=[], m=m, n=n, k=k, bs=bs, L=H0.L)
for alpha in np.arange(0, 26, 4.0):
    H = H0.scale_columns(10.0 ** (-alpha * u))
    A = H.dense()
    s = np.linalg.svd(A, compute_uv=False)
    nA = s[0]
    b = rng.standard_normal(m)
    xq = minnorm_qr(A, b)
    bw = lambda x: float(np.linalg.norm(b - A @ x) / (nA * np.linalg.norm(x) + np.linalg.norm(b)))
    r = dict(alpha=float(alpha), kappa=float(s[0] / s[-1]), qr_bwd=bw(xq))
    for key, root in [('ulv', 'factored'), ('ulv_pinv', 'pinv')]:
        F = MinNormFactor(H, slack='relax', root=root)
        x = F.solve(b)
        r[key] = dict(fwd=rel(x, xq), bwd=bw(x))
        r['root_shape'] = list(F.root_shape)
        R = F.root
        r['root_kappa'] = float(np.linalg.cond(R))
    out['rows'].append(r)
    print('alpha=%4.1f kappa=%.2e root %s cond %.1e | ulv fwd %.1e bwd %.1e | pinv-root fwd %.1e bwd %.1e | qr bwd %.1e'
          % (alpha, r['kappa'], r['root_shape'], r['root_kappa'], r['ulv']['fwd'], r['ulv']['bwd'],
             r['ulv_pinv']['fwd'], r['ulv_pinv']['bwd'], r['qr_bwd']), flush=True)
dump('m2b_deep', out)
