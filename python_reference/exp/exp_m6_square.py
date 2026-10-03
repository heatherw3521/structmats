"""E6 addendum: m/n = 1 (square) with the relaxed rule only.  (The MATLAB rule would
collapse every level -- slack l + n <= p fails at every square leaf -- and end in a dense
16384 x 16384 pinv; the repo routes square systems to hss_ulvvecsolve instead.)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from hssmn import *
from hssmn.refs import rel
from common import dump, load, timeit
rng = np.random.default_rng(55)
d = load('m5_sweeps'); d.pop('_machine', None); d.pop('_time', None)
n = 2 ** 14
H = random_hss(n, n, 16, 64, rng)
z = rng.standard_normal(n); xs = H.rmatvec(z); b = H.matvec(xs)
tf, _, F = timeit(lambda: MinNormFactor(H, slack='relax'), reps=3, warmup=0)
ts, _, x = timeit(lambda: F.solve(b), reps=3)
row = dict(ratio=1.0, m=n, n=n, L=H.L, slack_ok=bool(slack_ok(H)), relax=dict(t_factor=tf, t_solve=ts, err=rel(x, xs), res=rel(H.matvec(x), b)),
           collapse=None, note='MATLAB rule would collapse all levels (dense 16384x16384); not run')
d['aspect'] = [r for r in d['aspect'] if r['ratio'] != 1.0] + [row]
print(row)
dump('m5_sweeps', d)
