"""A2: 1-D deconvolution (decimated Gaussian blur), exact and noisy data.
A3: general (non-smooth) Toeplitz via the Cauchy-like transform: ranks, accuracy, scaling.
A4: 2-D blur (higher dimensions): HSS rank growth with N, lexicographic vs Morton order.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from hssmn import *
from hssmn.families import *
from hssmn.refs import rel, minnorm_qr
from hssmn.hss import matlab_tree
from common import dump, timeit

rng = np.random.default_rng(21)
which = sys.argv[1:] or ['a2', 'a3', 'a4']
out = {}


def test_signal(n):
    t = np.linspace(0, 1, n)
    s = (np.exp(-((t - 0.2) / 0.03) ** 2) + 0.6 * (np.abs(t - 0.5) < 0.08)
         + 0.8 * np.exp(-((t - 0.78) / 0.01) ** 2) + 0.3 * np.sin(6 * np.pi * t) * (t > 0.85))
    return t, s


if 'a2' in which:
    n = 2048
    sp = conv_spec(n, stride=2, kind='gauss', sigma=3.0)
    L, rb, cb = matlab_tree(sp.m, sp.n, 64)
    H = hss_build(sp, L, rb, cb, 1e-14, mode='proxy')
    A = sp.dense()
    sv = np.linalg.svd(A, compute_uv=False)
    t, x0 = test_signal(n)
    b = conv_apply(sp, x0)
    F = MinNormFactor(H)
    x_mn = F.solve(b)
    x_ref = minnorm_qr(A, b)
    noise = 1e-3 * np.linalg.norm(b) / np.sqrt(sp.m) * rng.standard_normal(sp.m)
    x_noisy = F.solve(b + noise)
    out['a2'] = dict(m=sp.m, n=n, L=H.L, maxrank=H.max_rank(), kappa=float(sv[0] / sv[-1]), sv=sv,
                     err_vs_dense=rel(x_mn, x_ref), res=rel(A @ x_mn, b),
                     rel_dist_truth_exact=rel(x_mn, x0), rel_dist_truth_noisy=rel(x_noisy, x0),
                     noise_level=1e-3, t=t, x0=x0, x_mn=x_mn, x_noisy=x_noisy, b=b, rpos=sp.rpos.real)
    print('A2 m=%d n=%d kappa=%.2e maxrank=%d err vs dense %.1e | ||x_mn-x0||/||x0|| exact %.2f noisy %.2e'
          % (sp.m, n, sv[0] / sv[-1], H.max_rank(), out['a2']['err_vs_dense'], out['a2']['rel_dist_truth_exact'],
             out['a2']['rel_dist_truth_noisy']), flush=True)
    dump('a2_a3_a4_' + '_'.join(which), out)

if 'a3' in which:
    rows = []
    for e in range(9, 17):
        n = 2 ** e
        m = n // 2
        t = rng.standard_normal(m + n - 1) / np.sqrt(n)
        tc = t[m - 1::-1].copy(); tr = t[m - 1:].copy(); tc[0] = tr[0]
        t0 = time.perf_counter()
        TC = ToeplitzCauchy(tc, tr)
        sp = TC.spec()
        L, rb, cb = matlab_tree(m, n, 64)
        H = hss_build(sp, L, rb, cb, 1e-12, mode='proxy')
        tb = time.perf_counter() - t0
        b = rng.standard_normal(m)
        tf, _, F = timeit(lambda: MinNormFactor(H), reps=3 if e < 15 else 1)
        def full_solve():
            y = F.solve(TC.Fm(b.astype(complex)))
            return TC.FnD_inv(y)
        ts, _, x = timeit(full_solve, reps=3)
        row = dict(m=m, n=n, L=H.L, rank_C=H.max_rank(), t_build=tb, t_factor=tf, t_solve_incl_fft=ts)
        # manufactured check in C coordinates (exact for the HSS matrix, any size)
        z = rng.standard_normal(m) + 1j * rng.standard_normal(m)
        ys = H.rmatvec(z); bb = H.matvec(ys)
        row['err_manuf'] = rel(F.solve(bb), ys)
        if e <= 12:
            T = TC.T_dense().real
            xr = minnorm_qr(T, b)
            row['err_vs_dense_T'] = rel(x, xr)
            row['res_T'] = rel(T @ x, b)
            row['imag_frac'] = float(np.linalg.norm(x.imag) / np.linalg.norm(x))
            LT, rbT, cbT = matlab_tree(m, n, 64)
            HT = hss_build(KernelSpec(m, n, lambda I, J: T[np.ix_(I, J)], np.arange(m), np.arange(n)),
                           LT, rbT, cbT, 1e-12, mode='full')
            row['rank_T'] = HT.max_rank()
        rows.append(row)
        print('A3', row, flush=True)
    out['a3'] = rows
    dump('a2_a3_a4_' + '_'.join(which), out)

if 'a4' in which:
    def morton(n1):
        idx = np.zeros((n1, n1), dtype=np.int64)
        for i in range(n1):
            for j in range(n1):
                v, b = 0, 0
                ii, jj = i, j
                while ii or jj:
                    v |= ((jj & 1) << (2 * b)) | ((ii & 1) << (2 * b + 1))
                    ii >>= 1; jj >>= 1; b += 1
                idx[i, j] = v
        return idx
    rows = []
    sig, w = 1.5, 5
    for n1 in [16, 32, 64, 128]:
        n2 = n1 // 2
        ii, jj = np.meshgrid(np.arange(n1), np.arange(n1), indexing='ij')
        oi, oj = np.meshgrid(2 * np.arange(n2) + 0.5, 2 * np.arange(n2) + 0.5, indexing='ij')
        for order in ['lex', 'morton']:
            if order == 'lex':
                pc = np.argsort((ii * n1 + jj).ravel()); pr = np.argsort((oi * n2 + oj).ravel())
            else:
                pc = np.argsort(morton(n1).ravel()); pr = np.argsort(morton(n2).ravel())
            ci, cj = ii.ravel()[pc], jj.ravel()[pc]
            ri, rj = oi.ravel()[pr], oj.ravel()[pr]
            def ent(I, J):
                d2 = (ri[I][:, None] - ci[J][None, :]) ** 2 + (rj[I][:, None] - cj[J][None, :]) ** 2
                out_ = np.exp(-d2 / (2 * sig ** 2))
                out_[d2 > w ** 2] = 0.0
                return out_
            m, n = n2 * n2, n1 * n1
            spec = KernelSpec(m, n, ent, np.zeros(m), np.zeros(n))
            L, rb, cb = matlab_tree(m, n, 32)
            t0 = time.perf_counter()
            H = hss_build(spec, L, rb, cb, 1e-10, mode='full')
            tb = time.perf_counter() - t0
            ranks = [max(max(u.shape[1] for u in H.U[l]), max(v.shape[1] for v in H.V[l])) for l in range(1, H.L + 1)]
            z = rng.standard_normal(m); xs = H.rmatvec(z); b = H.matvec(xs)
            tf, _, F = timeit(lambda: MinNormFactor(H), reps=3 if n <= 4096 else 1)
            x = F.solve(b)
            row = dict(n1=n1, N=n, m=m, order=order, L=H.L, maxrank=H.max_rank(), ranks_by_level=ranks,
                       t_build=tb, t_factor=tf, err=rel(x, xs), storage=H.storage())
            rows.append(row)
            print('A4', row, flush=True)
    out['a4'] = rows
    dump('a2_a3_a4_' + '_'.join(which), out)
