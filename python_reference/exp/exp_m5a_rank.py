"""E5a (revised): factor time vs HSS rank k at two FIXED leaf sizes.
F2, n = 2^16, m = n/2.  blocksize 64 -> 64 x 128 leaves (k = 2..32),
blocksize 128 -> 128 x 256 leaves (k = 2..64).  Compared with the flop count of
memo 1, Sec. 5 (hssmn.cost.factor_flops), scaled once per leaf size."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from hssmn import *
from hssmn.refs import rel
from hssmn.cost import factor_flops
from common import dump, load, timeit

rng = np.random.default_rng(55)
n = 2 ** 16; m = n // 2
rows = []
for bs, ks in [(64, [2, 4, 8, 16, 24, 32]), (128, [2, 4, 8, 16, 32, 48, 64])]:
    for k in ks:
        H = random_hss(m, n, k, bs, rng)
        z = rng.standard_normal(m); xs = H.rmatvec(z); b = H.matvec(xs)
        tf, _, F = timeit(lambda: MinNormFactor(H), reps=3, warmup=1)
        ts, _, x = timeit(lambda: F.solve(b), reps=3, warmup=1)
        n0 = m // 2 ** H.L
        rows.append(dict(k=k, bs=bs, L=H.L, n0=n0, t_factor=tf, t_solve=ts, err=rel(x, xs),
                         flops=factor_flops(m, n, H.L, k)))
        print('bs=%3d n0=%3d k=%2d factor %.3fs solve %.4fs flops %.2e err %.1e' % (bs, n0, k, tf, ts, rows[-1]['flops'], rows[-1]['err']), flush=True)
        del H, F
d = load('m5_sweeps')
d.pop('_machine', None); d.pop('_time', None)
d['rank_old'] = d.get('rank_old', d['rank'])
d['rank'] = rows
dump('m5_sweeps', d)
