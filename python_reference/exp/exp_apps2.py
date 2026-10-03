"""Memo 2 applications.

BP   sparse-spike deconvolution (Ricker wavelet, decimated), basis pursuit
     min ||x||_1 s.t. Hx = b by (i) IRLS: inner loop = diagonal-weight min-norm,
     (ii) Douglas-Rachford: inner step = min-norm projection, H factored ONCE.
DB   Tikhonov deblurring with noise: lambda sweep, discrepancy principle, L-curve.
MRI  k-space (Cartesian-grid) diagonal weights for non-Cartesian 1-D MRI:
     unweighted vs prior-weighted min-norm and weighted Tikhonov.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from hssmn import *
from hssmn.families import *
from hssmn.refs import rel, minnorm_qr
from hssmn.hss import matlab_tree, tree_from_leaves
from common import dump

rng = np.random.default_rng(41)
which = sys.argv[1:] or ['bp', 'db', 'mri']
out = {}
tag = '_'.join(which)

if 'bp' in which:
    n = 4096
    sp = conv_spec(n, stride=2, kind='ricker', sigma=2.5)
    L, rb, cb = matlab_tree(sp.m, sp.n, 64)
    H = hss_build(sp, L, rb, cb, 1e-14, mode='proxy')
    m = sp.m
    # sparse spike train, min separation 24 samples, away from the boundary
    s = 60
    pos = []
    while len(pos) < s:
        p = int(rng.integers(3 * sp.w, n - 3 * sp.w))
        if all(abs(p - q) >= 24 for q in pos):
            pos.append(p)
    x0 = np.zeros(n)
    x0[pos] = rng.choice([-1, 1], s) * (0.5 + rng.random(s))
    b = conv_apply(sp, x0)
    hist = dict(irls=[], dr=[])
    # minimum-norm solution for reference
    F0 = MinNormFactor(H)
    x_mn = F0.solve(b)
    # (i) IRLS (Daubechies-DeVore-Fornasier-Guntuerk 2010, eps-decreasing)
    x = x_mn.copy()
    eps = 1.0
    t0 = time.perf_counter()
    for it in range(60):
        w = (x ** 2 + eps ** 2) ** (-0.25)        # L = diag(w): sum w_j^2 x_j^2 = sum x_j^2/sqrt(x_j^2+eps^2)
        x = weighted_minnorm(H, b, w=w)
        r_sorted = np.sort(np.abs(x))[::-1]
        eps = min(eps, r_sorted[s] / n) if r_sorted[s] > 0 else eps / 10
        eps = max(eps, 1e-10)
        hist['irls'].append(dict(it=it + 1, t=time.perf_counter() - t0, err=rel(x, x0),
                                 res=rel(H.matvec(x), b), l1=float(np.sum(np.abs(x))), eps=float(eps)))
        if rel(x, x0) < 1e-10:
            break
    x_irls = x
    # (ii) Douglas-Rachford on  ||x||_1 + indicator{Hx=b}, projection P(v) = v + H^+(b - Hv)
    gamma = 0.05
    t0 = time.perf_counter()
    F = MinNormFactor(H)
    t_fac = time.perf_counter() - t0
    proj = lambda v: v + F.solve(b - H.matvec(v))
    soft = lambda v, g: np.sign(v) * np.maximum(np.abs(v) - g, 0)
    zz = x_mn.copy()
    for it in range(3000):
        xp = proj(zz)
        zz = zz + soft(2 * xp - zz, gamma) - xp
        if (it + 1) % 10 == 0 or it < 10:
            hist['dr'].append(dict(it=it + 1, t=t_fac + time.perf_counter() - t0 - t_fac, err=rel(xp, x0),
                                   res=rel(H.matvec(xp), b), l1=float(np.sum(np.abs(xp)))))
        if rel(xp, x0) < 1e-8:
            hist['dr'].append(dict(it=it + 1, t=time.perf_counter() - t0, err=rel(xp, x0),
                                   res=rel(H.matvec(xp), b), l1=float(np.sum(np.abs(xp)))))
            break
    x_dr = xp
    # timing of one factorization vs one solve
    t0 = time.perf_counter(); Ft = MinNormFactor(H); t_factor = time.perf_counter() - t0
    t0 = time.perf_counter()
    for _ in range(20):
        Ft.solve(b)
    t_solve = (time.perf_counter() - t0) / 20
    out['bp'] = dict(m=m, n=n, s=s, L=H.L, maxrank=H.max_rank(), x0=x0, x_mn=x_mn, x_irls=x_irls, x_dr=x_dr,
                     hist=hist, t_factor=t_factor, t_solve=t_solve, gamma=gamma,
                     err_mn=rel(x_mn, x0), err_irls=rel(x_irls, x0), err_dr=rel(x_dr, x0))
    print('BP m=%d n=%d s=%d | min-norm err %.2e | IRLS %d its err %.1e | DR %d its err %.1e | factor %.3fs solve %.4fs'
          % (m, n, s, rel(x_mn, x0), len(hist['irls']), rel(x_irls, x0), hist['dr'][-1]['it'], rel(x_dr, x0),
             t_factor, t_solve), flush=True)
    dump('apps2_' + tag, out)

if 'db' in which:
    n = 4096
    sp = conv_spec(n, stride=2, kind='gauss', sigma=3.0)
    L, rb, cb = matlab_tree(sp.m, sp.n, 64)
    H = hss_build(sp, L, rb, cb, 1e-14, mode='proxy')
    t = np.linspace(0, 1, n)
    x0 = (np.exp(-((t - 0.2) / 0.03) ** 2) + 0.6 * (np.abs(t - 0.5) < 0.08)
          + 0.8 * np.exp(-((t - 0.78) / 0.01) ** 2) + 0.3 * np.sin(6 * np.pi * t) * (t > 0.85))
    b0 = conv_apply(sp, x0)
    noise_rel = 1e-2
    e = noise_rel * np.linalg.norm(b0) / np.sqrt(sp.m) * rng.standard_normal(sp.m)
    b = b0 + e
    lams = np.logspace(-5, 0.5, 23)
    rows = []
    t0 = time.perf_counter()
    for lam in lams:
        x = tikhonov(H, b, lam)
        rows.append(dict(lam=float(lam), err=rel(x, x0), res=float(np.linalg.norm(H.matvec(x) - b)),
                         xnorm=float(np.linalg.norm(x))))
    t_sweep = time.perf_counter() - t0
    nrm_e = float(np.linalg.norm(e))
    # discrepancy principle: largest lambda with residual <= 1.01*||e||
    ok = [r for r in rows if r['res'] <= 1.01 * nrm_e]
    lam_dp = max(r['lam'] for r in ok) if ok else None
    best = min(rows, key=lambda r: r['err'])
    x_dp = tikhonov(H, b, lam_dp) if lam_dp else None
    x_best = tikhonov(H, b, best['lam'])
    x_mn = MinNormFactor(H).solve(b)
    out['db'] = dict(m=sp.m, n=n, noise_rel=noise_rel, noise_norm=nrm_e, rows=rows, lam_dp=lam_dp,
                     lam_best=best['lam'], err_best=best['err'], err_dp=rel(x_dp, x0) if lam_dp else None,
                     err_mn=rel(x_mn, x0), t_per_lambda=t_sweep / len(lams), t=t, x0=x0, x_dp=x_dp, x_best=x_best,
                     x_mn=x_mn)
    print('DB m=%d n=%d | min-norm err %.2e | best lam %.1e err %.3f | discrepancy lam %s err %s | %.2fs per lambda'
          % (sp.m, n, rel(x_mn, x0), best['lam'], best['err'], lam_dp, out['db']['err_dp'], t_sweep / len(lams)), flush=True)
    dump('apps2_' + tag, out)

if 'mri' in which:
    N = 2048
    tgrid = np.arange(-N // 2, N // 2)
    # 1-D phantom (piecewise constant + smooth), centred image grid
    x0 = (1.0 * (np.abs(tgrid) < 600) + 0.5 * (np.abs(tgrid - 150) < 120) - 0.4 * (np.abs(tgrid + 260) < 60)
          + 0.3 * np.exp(-((tgrid - 420) / 25.0) ** 2))
    # variable-density non-Cartesian k-space samples xi in [-1/2,1/2): denser at low |k|
    m = int(0.45 * N)
    u = np.sort(rng.random(m))
    # inverse-CDF of density ~ 1/(1+|xi|/0.12) on [-1/2,1/2)
    a = 0.12
    Fcdf = lambda r: np.sign(r) * a * np.log1p(np.abs(r) / a)
    tot = Fcdf(0.5)
    q = (u - 0.5) * 2 * tot
    xi = np.sign(q) * a * np.expm1(np.abs(q) / a)
    xi = np.clip(xi, -0.5 + 1e-9, 0.5 - 1e-9)
    xi = xi + 1e-6 * rng.standard_normal(m) / N
    # F5 form: A_jt = e^{-2 pi i xi_j t} = e^{2 pi i t x'_j}, x' = -xi mod 1
    xp = np.sort(np.mod(-xi, 1.0))
    sp = nudft_spec(xp, N)
    A = nudft_vandermonde(xp, N)          # columns t = -N/2..N/2-1 (image grid)
    d = A @ x0
    Fu = unitary_dft(N)                   # y = F x : Cartesian k-space (grid l/N)
    y0 = Fu @ x0
    # geometry-aligned tree (rows split where the k-space grid splits)
    L, rb, cb = matlab_tree(m, N, 32)
    gridpos = np.arange(N) / N
    rbL = np.searchsorted(xp, gridpos[cb[L][:-1]]); rbL = np.append(rbL, m); rbL[0] = 0
    Lg, rbg, cbg = tree_from_leaves(rbL, cb[L])
    H = hss_build(sp, Lg, rbg, cbg, 1e-12, mode='proxy')
    Hidx = hss_build(sp, L, rb, cb, 1e-12, mode='proxy')
    # Cartesian frequency of column l: f_l = l for l < N/2, l - N otherwise  (grid l/N, periodic)
    lidx = np.arange(N)
    kfreq = np.where(lidx < N // 2, lidx, lidx - N)
    # prior power spectrum of the phantom class ~ (1 + |k|/k0)^{-2}
    k0 = 8.0
    wts = (1 + np.abs(kfreq) / k0)       # L = diag(w), w = P^{-1/2}
    y_mn = MinNormFactor(H).solve(d)
    y_w = weighted_minnorm(H, d, w=wts)
    k0sweep = []
    for kk0 in [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 1e9]:
        yk = weighted_minnorm(H, d, w=(1 + np.abs(kfreq) / kk0))
        k0sweep.append(dict(k0=kk0, err=rel(Fu.conj().T @ yk, x0)))
    x_mn = Fu.conj().T @ y_mn
    x_w = Fu.conj().T @ y_w
    # noisy data: weighted Tikhonov  min ||Cy - d_n||^2 + lam^2 ||W y||^2
    sig = 2e-2 * np.linalg.norm(d) / np.sqrt(m)
    dn = d + sig * (rng.standard_normal(m) + 1j * rng.standard_normal(m)) / np.sqrt(2)
    tik = []
    for lam in np.logspace(-4, 0, 13):
        yt = tikhonov(H, dn, lam, w=wts)
        tik.append(dict(lam=float(lam), err=rel(Fu.conj().T @ yt, x0)))
    bestt = min(tik, key=lambda r: r['err'])
    y_t = tikhonov(H, dn, bestt['lam'], w=wts)
    y_mn_noisy = MinNormFactor(H).solve(dn)
    C = sp.dense()
    out['mri'] = dict(m=m, N=N, rank_geo=H.max_rank(), rank_idx=Hidx.max_rank(),
                      err_minnorm=rel(x_mn, x0), err_weighted=rel(x_w, x0),
                      err_minnorm_noisy=rel(Fu.conj().T @ y_mn_noisy, x0), err_weighted_tik=bestt['err'],
                      lam_best=bestt['lam'], tik=tik, k0sweep=k0sweep, k0=k0,
                      res_w=rel(C @ y_w, d), check_y0=rel(C @ y0, d),
                      t=tgrid, x0=x0, x_mn=x_mn.real, x_w=x_w.real, x_tik=(Fu.conj().T @ y_t).real,
                      xi=xi)
    print('MRI m=%d N=%d rank geo %d idx %d | image err: min-norm %.3f weighted %.3f | noisy: min-norm %.3f weighted-Tikhonov %.3f (lam %.1e)'
          % (m, N, H.max_rank(), Hidx.max_rank(), out['mri']['err_minnorm'], out['mri']['err_weighted'],
             out['mri']['err_minnorm_noisy'], bestt['err'], bestt['lam']), flush=True)
    dump('apps2_' + tag, out)
