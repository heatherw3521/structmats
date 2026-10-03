"""E3/E4: large-scale accuracy (manufactured min-norm solutions, no dense
matrix anywhere) and timing/scaling.

Manufactured solution:  z ~ N(0,I_m),  x* = H^* z,  b = H x*.  Since x* lies
in Range(H^*) and solves Hx = b, it IS the minimum-norm solution of the HSS
system, so ||x - x*||/||x*|| measures the solver alone, at any size, using
only O(N) HSS products.

Run with OPENBLAS_NUM_THREADS=1 (timings are single-threaded, median of reps).
"""
import sys, os, time, gc
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from hssmn import *
from hssmn.families import *
from hssmn.refs import rel, minnorm_qr, minnorm_cgne
from hssmn.hss import matlab_tree
from common import dump, timeit, load

rng = np.random.default_rng(3)
QUICK = '--quick' in sys.argv
maxexp = 14 if QUICK else 20


def build(fam, n):
    """returns (H, info) -- never forms a dense matrix"""
    if fam == 'F1':
        m = n // 2
        return lrbd_hss(m, n, 3, 16, rng), {}
    if fam == 'F2':
        m = n // 2
        return random_hss(m, n, 10, 64, rng), {}
    if fam == 'F3':
        M = n // 3 + 1
        sp = cauchy_spec(M, 3)
        L, rb, cb = matlab_tree(sp.m, sp.n, 64)
        t0 = time.perf_counter()
        H = hss_build(sp, L, rb, cb, 1e-9, mode='proxy')
        return H, dict(build=time.perf_counter() - t0)
    if fam == 'F4':
        sp = conv_spec(n, stride=2, kind='gauss', sigma=2.0)
        L, rb, cb = matlab_tree(sp.m, sp.n, 64)
        t0 = time.perf_counter()
        H = hss_build(sp, L, rb, cb, 1e-13, mode='proxy')
        return H, dict(build=time.perf_counter() - t0)
    if fam == 'F5':
        m = (3 * n) // 4
        sp = nudft_spec(nudft_points(m, rng), n)
        L, rb, cb = matlab_tree(sp.m, sp.n, 64)
        t0 = time.perf_counter()
        H = hss_build(sp, L, rb, cb, 1e-10, mode='proxy')
        return H, dict(build=time.perf_counter() - t0)
    raise ValueError(fam)


ranges = {'F1': range(10, maxexp + 1), 'F2': range(10, maxexp + 1), 'F3': range(10, min(maxexp, 19) + 1),
          'F4': range(10, min(maxexp, 19) + 1), 'F5': range(10, min(maxexp, 18) + 1)}
only = [a for a in sys.argv[1:] if a.startswith('F')]
resname = 'm3_scaling' + ('_quick' if QUICK else '')
try:
    out = load(resname) if only else {}
except FileNotFoundError:
    out = {}
out.pop('_machine', None); out.pop('_time', None)
for fam, exps in ranges.items():
    if only and fam not in only:
        continue
    out[fam] = []
    for e in exps:
        n = 2 ** e
        H, info = build(fam, n)
        m, n = H.shape
        cplx = np.iscomplexobj(np.zeros(0, dtype=H.dtype))
        z = rng.standard_normal(m) + (1j * rng.standard_normal(m) if cplx else 0)
        xs = H.rmatvec(z)
        b = H.matvec(xs)
        reps = 5 if n <= 2 ** 16 else 3
        tf, tfs, F = timeit(lambda: MinNormFactor(H, slack='relax'), reps=reps, warmup=1 if n <= 2 ** 16 else 0)
        ts, tss, x = timeit(lambda: F.solve(b), reps=reps, warmup=1)
        tm, _, _ = timeit(lambda: H.matvec(xs), reps=reps, warmup=1)
        r = dict(m=m, n=n, L=H.L, maxrank=H.max_rank(), storage=H.storage(),
                 t_factor=tf, t_factor_all=tfs, t_solve=ts, t_matvec=tm,
                 res=rel(H.matvec(x), b), err=rel(x, xs), **info)
        out[fam].append(r)
        print('%s n=%8d m=%8d L=%2d r=%3d factor %.3fs solve %.4fs matvec %.4fs | res %.1e err %.1e %s' %
              (fam, n, m, H.L, r['maxrank'], tf, ts, tm, r['res'], r['err'],
               ('build %.2fs' % info['build']) if 'build' in info else ''), flush=True)
        del H, F, x
        gc.collect()
    dump(resname, out)

# ---- dense comparison (QR-based min-norm) and CGNE, moderate sizes
dense = []
for e in range(9, 14):
    n = 2 ** e
    H = random_hss(n // 2, n, 10, 64, rng)
    A = H.dense()
    b = rng.standard_normal(n // 2)
    reps = 3 if e < 13 else 1
    td, _, xd = timeit(lambda: minnorm_qr(A, b), reps=reps, warmup=1 if e < 13 else 0)
    tf, _, F = timeit(lambda: MinNormFactor(H), reps=3, warmup=1)
    ts, _, x = timeit(lambda: F.solve(b), reps=3, warmup=1)
    tc, _, (xc, its, info) = timeit(lambda: minnorm_cgne(H, b, rtol=1e-12), reps=1, warmup=0)
    dense.append(dict(n=n, t_dense=td, t_factor=tf, t_solve=ts, t_cgne=tc, cgne_its=its,
                      err_hss_vs_dense=rel(x, xd), err_cgne=rel(xc, xd)))
    print('dense n=%6d qr %.3fs | hss factor %.3fs solve %.4fs | cgne %.3fs (%d its)' % (n, td, tf, ts, tc, its), flush=True)
    del A
out['dense_F2'] = dense
dump(resname, out)
