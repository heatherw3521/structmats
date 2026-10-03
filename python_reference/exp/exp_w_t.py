"""Memo 2 accuracy/timing experiments.

W1  weighted min-norm (diagonal L = diag(w), w_j = 10^{alpha u_j}, u ~ U[-1/2,1/2];
    leaf-block L_tau SPD) vs a dense backward-stable reference, families F2, F3, F4
W2  large-scale weighted, manufactured solution x* = L^{-1} L^{-*} H^* z
T1  Tikhonov lambda sweep for wide / square / tall H vs the dense SVD filter solution;
    general form (row weights S, column weights L)
T2  large-scale Tikhonov, manufactured  x* = H^* z,  b = H x* + lam^2 z
T3  timing: Tikhonov (augmented) vs plain min-norm vs N
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from hssmn import *
from hssmn.families import *
from hssmn.refs import rel, minnorm_qr, tikhonov_svd
from hssmn.hss import matlab_tree
from common import dump, timeit

rng = np.random.default_rng(31)
which = sys.argv[1:] or ['w1', 'w2', 't1', 't2', 't3']
out = {}
tag = '_'.join(which)


def proto(n):
    """memo 1 timing protocol (E4): warm-up run, then median of 5 (n <= 2^16) or 3."""
    return (5, 1) if n <= 2 ** 16 else (3, 0)


def fam(name, size='small'):
    if name == 'F2':
        return random_hss(1024, 2048, 10, 32, rng)
    if name == 'F3':
        sp = cauchy_spec(683, 3)
        L, rb, cb = matlab_tree(sp.m, sp.n, 32)
        return hss_build(sp, L, rb, cb, 1e-12, mode='proxy')
    if name == 'F4':
        sp = conv_spec(2048, stride=2, kind='gauss', sigma=2.0)
        L, rb, cb = matlab_tree(sp.m, sp.n, 32)
        return hss_build(sp, L, rb, cb, 1e-14, mode='proxy')
    if name == 'F2sq':
        return random_hss(1536, 1536, 10, 32, rng)
    if name == 'F2tall':
        return random_hss(2048, 1024, 10, 32, rng)
    if name == 'F4sq':
        sp = conv_spec(1536, stride=1, kind='gauss', sigma=2.0, mode='same')
        L, rb, cb = matlab_tree(sp.m, sp.n, 32)
        return hss_build(sp, L, rb, cb, 1e-14, mode='proxy')
    raise ValueError(name)


if 'w1' in which:
    rows = []
    for name in ['F2', 'F3', 'F4']:
        H = fam(name)
        A = H.dense()
        m, n = A.shape
        b = rng.standard_normal(m)
        u = rng.random(n) - 0.5
        for alpha in [0, 2, 4, 6, 8, 10, 12]:
            w = 10.0 ** (alpha * u)
            x = weighted_minnorm(H, b, w=w)
            Aw = A / w[None, :]
            sv = np.linalg.svd(Aw, compute_uv=False)
            xr = minnorm_qr(Aw, b) / w
            rows.append(dict(family=name, kind='diag', alpha=alpha, m=m, n=n, kappa_HLinv=float(sv[0] / sv[-1]),
                             err=rel(x, xr), res=rel(A @ x, b),
                             wnorm_ratio=float(np.linalg.norm(w * x) / np.linalg.norm(w * xr))))
            print('W1', rows[-1], flush=True)
        # leaf-block weights: L_tau = chol(I + beta * K_tau)^T with a smoothing kernel K_tau
        L_ = H.L; cbL = H.cb[L_]
        for beta in [1.0, 1e2, 1e4]:
            Lb = []
            for i in range(2 ** L_):
                p = cbL[i + 1] - cbL[i]
                tt = np.arange(p)
                K = np.exp(-np.subtract.outer(tt, tt) ** 2 / (2 * 3.0 ** 2))
                Lb.append(np.linalg.cholesky(np.eye(p) + beta * K).T)
            x = weighted_minnorm(H, b, Lblocks=Lb)
            Lfull = np.zeros((n, n))
            for i in range(2 ** L_):
                Lfull[cbL[i]:cbL[i + 1], cbL[i]:cbL[i + 1]] = Lb[i]
            Li = np.linalg.inv(Lfull)
            ALi = A @ Li
            sv = np.linalg.svd(ALi, compute_uv=False)
            xr = Li @ minnorm_qr(ALi, b)
            rows.append(dict(family=name, kind='block', beta=beta, m=m, n=n, kappa_HLinv=float(sv[0] / sv[-1]),
                             err=rel(x, xr), res=rel(A @ x, b)))
            print('W1', rows[-1], flush=True)
    out['w1'] = rows
    dump('wt_' + tag, out)

if 'w2' in which:
    rows = []
    for e in range(12, 20):
        n = 2 ** e
        H = random_hss(n // 2, n, 10, 64, rng)
        w = 10.0 ** (6 * (rng.random(n) - 0.5))
        z = rng.standard_normal(n // 2)
        xs = H.scale_columns(1 / w).rmatvec(z) / w      # L^{-1} (H L^{-1})^* z
        b = H.matvec(xs)
        reps, wu = proto(n)
        tw, _, x = timeit(lambda: weighted_minnorm(H, b, w=w), reps=reps, warmup=wu)
        tu, _, _ = timeit(lambda: MinNormFactor(H).solve(b), reps=reps, warmup=wu)
        rows.append(dict(n=n, t=tw, t_unweighted=tu, err=rel(x, xs), res=rel(H.matvec(x), b)))
        print('W2', rows[-1], flush=True)
    out['w2'] = rows
    dump('wt_' + tag, out)

if 't1' in which:
    rows = []
    for name in ['F2', 'F2sq', 'F2tall', 'F4', 'F4sq']:
        H = fam(name)
        A = H.dense()
        m, n = A.shape
        U, s, Vh = np.linalg.svd(A, full_matrices=False)
        b = rng.standard_normal(m)
        for lr in np.arange(-8, 3, 1.0):
            lam = 10.0 ** lr * s[0]
            x = tikhonov(H, b, lam)
            xr = Vh.conj().T @ ((s / (s ** 2 + lam ** 2)) * (U.conj().T @ b))
            # residual of the regularised normal equations, relative
            g = A.conj().T @ (A @ x - b) + lam ** 2 * x
            rows.append(dict(family=name, m=m, n=n, lam_rel=float(10.0 ** lr), err=rel(x, xr),
                             kappa_A=float(s[0] / s[-1]), kappa_aug=float(np.sqrt(s[0] ** 2 + lam ** 2) / np.sqrt(s[-1] ** 2 + lam ** 2)) if m <= n else float(np.sqrt(s[0] ** 2 + lam ** 2) / lam),
                             normal_eq_res=float(np.linalg.norm(g) / (np.linalg.norm(A.conj().T @ b)))))
            print('T1', name, '%dx%d lam/||A||=%.0e err %.1e' % (m, n, 10.0 ** lr, rows[-1]['err']), flush=True)
    # general form
    H = fam('F2')
    A = H.dense(); m, n = A.shape
    b = rng.standard_normal(m)
    for lr in [-4, -2, 0]:
        lam = 10.0 ** lr
        w = 10.0 ** (2 * (rng.random(n) - 0.5)); sr = 10.0 ** (2 * (rng.random(m) - 0.5))
        x = tikhonov(H, b, lam, w=w, s=sr)
        # dense: stacked least squares [S A; lam W] x ~ [S b; 0]
        M = np.vstack([sr[:, None] * A, lam * np.diag(w)])
        rhs = np.concatenate([sr * b, np.zeros(n)])
        xr = np.linalg.lstsq(M, rhs, rcond=None)[0]
        rows.append(dict(family='F2 general form', m=m, n=n, lam_rel=lam, err=rel(x, xr)))
        print('T1 general', rows[-1], flush=True)
    out['t1'] = rows
    dump('wt_' + tag, out)

if 't2' in which or 't3' in which:
    rows = []
    for e in range(10, 20):
        n = 2 ** e
        H = random_hss(n // 2, n, 10, 64, rng)
        lam = 1e-2
        z = rng.standard_normal(n // 2)
        xs = H.rmatvec(z)
        b = H.matvec(xs) + lam ** 2 * z
        reps, wu = proto(n)
        tt, _, x = timeit(lambda: tikhonov(H, b, lam), reps=reps, warmup=wu)
        tm, _, xm = timeit(lambda: MinNormFactor(H).solve(b), reps=reps, warmup=wu)
        rows.append(dict(n=n, lam=lam, t_tikhonov=tt, t_minnorm=tm, err=rel(x, xs)))
        print('T2/T3', rows[-1], flush=True)
    # tall: Tikhonov is the only route (min-norm solver does not apply)
    for e in range(10, 19):
        n = 2 ** e
        H = random_hss(2 * n, n, 10, 64, rng)
        lam = 1e-2
        z = rng.standard_normal(2 * n)
        xs = H.rmatvec(z)
        b = H.matvec(xs) + lam ** 2 * z
        reps, wu = proto(n)
        tt, _, x = timeit(lambda: tikhonov(H, b, lam), reps=reps, warmup=wu)
        rows.append(dict(n=n, m=2 * n, shape='tall', lam=lam, t_tikhonov=tt, err=rel(x, xs)))
        print('T2 tall', rows[-1], flush=True)
    out['t2'] = rows
    dump('wt_' + tag, out)
