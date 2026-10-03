"""E5/E6: parameter sweeps (single-threaded timings, median of reps).

  E5a  time vs HSS rank k (F2, n = 2^16, m = n/2, blocksize 64 or 2k)
  E5b  time vs blocksize (F2, n = 2^16, k = 10)
  E6   aspect ratio m/n sweep (F2, n = 2^14, k = 16, blocksize 64):
       MATLAB slack handling ('collapse': merge the leaf level while some leaf
       violates l + n_tau <= p_tau) versus the relaxed rule p' = min(p, l + n)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from hssmn import *
from hssmn.refs import rel
from hssmn.hss import matlab_tree
from common import dump, timeit

rng = np.random.default_rng(5)
out = dict(rank=[], blocksize=[], aspect=[])

n = 2 ** 16
for k in [2, 4, 8, 16, 32, 48, 64]:
    bs = max(64, 2 * k)
    H = random_hss(n // 2, n, k, bs, rng)
    z = rng.standard_normal(n // 2); xs = H.rmatvec(z); b = H.matvec(xs)
    tf, _, F = timeit(lambda: MinNormFactor(H), reps=3)
    ts, _, x = timeit(lambda: F.solve(b), reps=3)
    out['rank'].append(dict(k=k, bs=bs, L=H.L, t_factor=tf, t_solve=ts, err=rel(x, xs), storage=H.storage()))
    print('k=%3d bs=%3d L=%d factor %.3fs solve %.4fs err %.1e' % (k, bs, H.L, tf, ts, rel(x, xs)), flush=True)

for bs in [16, 32, 64, 128, 256, 512]:
    H = random_hss(n // 2, n, 10, bs, rng)
    z = rng.standard_normal(n // 2); xs = H.rmatvec(z); b = H.matvec(xs)
    tf, _, F = timeit(lambda: MinNormFactor(H), reps=3)
    ts, _, x = timeit(lambda: F.solve(b), reps=3)
    out['blocksize'].append(dict(bs=bs, L=H.L, leaf=H.summary()['leaf_rows'] + H.summary()['leaf_cols'],
                                 t_factor=tf, t_solve=ts, err=rel(x, xs)))
    print('bs=%3d L=%2d factor %.3fs solve %.4fs err %.1e' % (bs, H.L, tf, ts, rel(x, xs)), flush=True)

n = 2 ** 14
for ratio in [0.1, 0.25, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]:
    m = int(round(ratio * n))
    H = random_hss(m, n, 16, 64, rng)
    z = rng.standard_normal(m); xs = H.rmatvec(z); b = H.matvec(xs)
    row = dict(ratio=ratio, m=m, n=n, L=H.L, slack_ok=bool(slack_ok(H)),
               leaf=H.summary()['leaf_rows'] + H.summary()['leaf_cols'])
    for slack in ['relax', 'collapse']:
        reps = 3
        tf, _, F = timeit(lambda: MinNormFactor(H, slack=slack), reps=reps, warmup=0)
        ts, _, x = timeit(lambda: F.solve(b), reps=3)
        row[slack] = dict(t_factor=tf, t_solve=ts, err=rel(x, xs), res=rel(H.matvec(x), b),
                          collapses=int(sum(F.collapses)), collapses_by_level=list(map(int, F.collapses)),
                          root=list(map(int, F.root_shape)), nlevels_processed=len(F.levels))
    out['aspect'].append(row)
    print('m/n=%.2f slack_ok=%s | relax %.3fs err %.1e | collapse %.3fs err %.1e collapses %s root %s' %
          (ratio, row['slack_ok'], row['relax']['t_factor'], row['relax']['err'], row['collapse']['t_factor'],
           row['collapse']['err'], row['collapse']['collapses_by_level'], row['collapse']['root']), flush=True)
    dump('m5_sweeps', out)
dump('m5_sweeps', out)
