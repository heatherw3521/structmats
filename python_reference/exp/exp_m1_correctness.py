"""E1: correctness of the min-norm solver against dense references, all families.

For each case:  kappa(H), residual, error vs the dense backward-stable min-norm
solution of the HSS matrix H itself (isolates solver error), error vs pinv of
the original matrix A (adds compression error), null-space fraction of x,
and the MATLAB-behaviour ('collapse') vs relaxed slack handling.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from hssmn import *
from hssmn.families import *
from hssmn.refs import *
from hssmn.hss import matlab_tree
from common import dump

rng = np.random.default_rng(2026)
rows = []


def run_case(name, family, H, A=None, transform=None, note=''):
    m, n = H.shape
    Ah = H.dense()
    cplx = np.iscomplexobj(Ah)
    b = rng.standard_normal(m) + (1j * rng.standard_normal(m) if cplx else 0)
    s = np.linalg.svd(Ah, compute_uv=False)
    kap = float(s[0] / s[-1])
    xq = minnorm_qr(Ah, b)
    out = dict(name=name, family=family, m=m, n=n, L=H.L, maxrank=H.max_rank(),
               leaf=H.summary()['leaf_rows'] + H.summary()['leaf_cols'], kappa=kap, note=note)
    for slack in ['collapse', 'relax']:
        t0 = time.perf_counter()
        F = MinNormFactor(H, slack=slack)
        x = F.solve(b)
        t1 = time.perf_counter() - t0
        out[slack] = dict(collapses=int(sum(F.collapses)), root=list(F.root_shape),
                          res=rel(Ah @ x, b), err_qr=rel(x, xq), null=nullspace_fraction(Ah, x),
                          normratio=float(np.linalg.norm(x) / np.linalg.norm(xq)), time=t1)
    if A is not None:
        out['comprerr'] = float(np.linalg.norm(Ah - A) / np.linalg.norm(A))
        xa = minnorm_qr(A, b)
        x = MinNormFactor(H, slack='relax').solve(b)
        out['err_vs_A'] = rel(x, xa)
    if transform is not None:
        out.update(transform(H, b))
    rows.append(out)
    r = out['relax']
    print('%-28s %5dx%-6d L=%d r=%2d kappa=%.1e | res %.1e err(QR of H) %.1e null %.1e | collapses(MATLAB)=%d root %s'
          % (name, m, n, H.L, H.max_rank(), kap, r['res'], r['err_qr'], r['null'],
             out['collapse']['collapses'], out['collapse']['root']), flush=True)


# ---- F1: LR + block diagonal (showcase rectsampler), exact HSS
for (m, n, k, bs) in [(40, 60, 6, 15), (320, 480, 3, 15), (2048, 4096, 3, 16)]:
    run_case('F1 lrbd k=%d bs=%d' % (k, bs), 'F1', lrbd_hss(m, n, k, bs, rng))
# ---- F2: random exact HSS (real, complex, aspect 0.9 -> slack failure)
run_case('F2 random k=10', 'F2', random_hss(1024, 2048, 10, 32, rng))
run_case('F2 random complex k=10', 'F2', random_hss(1024, 2048, 10, 32, rng, cplx=True))
run_case('F2 random m/n=0.9 k=8', 'F2', random_hss(1800, 2000, 8, 32, rng), note='slack fails')
# ---- F3: interlaced Cauchy (constructor rule 'full' and proxy build)
spec = cauchy_spec(1024, 3)
A = spec.dense()
L, rb, cb = matlab_tree(spec.m, spec.n, 32)
run_case('F3 cauchy (full ID)', 'F3', hss_build(spec, L, rb, cb, 1e-10, mode='full'), A)
run_case('F3 cauchy (proxy)', 'F3', hss_build(spec, L, rb, cb, 1e-10, mode='proxy'), A)
# ---- F4: decimated convolution
for kind, sig, st in [('gauss', 2.0, 2), ('ricker', 3.0, 2), ('gauss', 3.0, 3)]:
    sp = conv_spec(4096, stride=st, kind=kind, sigma=sig)
    A = sp.dense()
    L, rb, cb = matlab_tree(sp.m, sp.n, 32)
    run_case('F4 conv %s s=%g stride=%d' % (kind, sig, st), 'F4', hss_build(sp, L, rb, cb, 1e-13, mode='proxy'), A)
# ---- F5: NUDFT Cauchy-like (band-limited interpolation)
xs = nudft_points(1536, rng)
sp = nudft_spec(xs, 2048)
A = sp.dense()
L, rb, cb = matlab_tree(sp.m, sp.n, 32)
run_case('F5 nudft cauchy-like', 'F5', hss_build(sp, L, rb, cb, 1e-12, mode='proxy'), A)
# ---- F6: random Toeplitz via Cauchy-like transform
m6, n6 = 1024, 2048
t = rng.standard_normal(m6 + n6 - 1) / np.sqrt(n6)
tc = t[m6 - 1::-1].copy(); tr = t[m6 - 1:].copy(); tc[0] = tr[0]
TC = ToeplitzCauchy(tc, tr)
T = TC.T_dense().real
sp = TC.spec()
L, rb, cb = matlab_tree(m6, n6, 32)
HC = hss_build(sp, L, rb, cb, 1e-12, mode='proxy')
LT, rbT, cbT = matlab_tree(m6, n6, 32)
HT = hss_build(KernelSpec(m6, n6, lambda I, J: T[np.ix_(I, J)], np.arange(m6), np.arange(n6)), LT, rbT, cbT, 1e-12, mode='full')

def toep_back(H, b):
    bb = rng.standard_normal(m6)
    y = MinNormFactor(H, slack='relax').solve(TC.Fm(bb.astype(complex)))
    x = TC.FnD_inv(y)
    xr = minnorm_qr(T, bb)
    return dict(toeplitz_res=rel(T @ x, bb), toeplitz_err=rel(x, xr), toeplitz_imag=float(np.linalg.norm(x.imag) / np.linalg.norm(x)),
                rank_T_itself=HT.max_rank())

run_case('F6 toeplitz->cauchy-like', 'F6', HC, TC.C_dense_fft(), transform=toep_back)
print('F6: max HSS rank of T itself = %d, of C = %d' % (HT.max_rank(), HC.max_rank()))
dump('m1_correctness', dict(rows=rows))
