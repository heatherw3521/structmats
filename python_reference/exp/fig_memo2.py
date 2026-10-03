"""Figures for memo 2 from results/*.json."""
import sys
import numpy as np
import matplotlib.pyplot as plt
from common import style, series, ref_slope, save, load, C, MK, INK2, MUTED

style()
U = np.finfo(float).eps / 2
which = sys.argv[1:] or ['w', 't', 'bp', 'db', 'mri']

if 'w' in which:
    d = load('wt_w1')['w1']
    fig, ax = plt.subplots(1, 1, figsize=(4.2, 3.1))
    names = {'F2': 'F2 random HSS', 'F3': 'F3 Cauchy', 'F4': 'F4 Gaussian conv.'}
    for i, f in enumerate(['F2', 'F3', 'F4']):
        rows = [r for r in d if r['family'] == f and r['kind'] == 'diag']
        k = np.array([r['kappa_HLinv'] for r in rows]); e = np.array([r['err'] for r in rows])
        o = np.argsort(k)
        series(ax, k[o], e[o], i, names[f] + ', diagonal $L$')
        rb = [r for r in d if r['family'] == f and r['kind'] == 'block']
        ax.plot([r['kappa_HLinv'] for r in rb], [r['err'] for r in rb], linestyle='none', marker=MK[i],
                markerfacecolor='white', markeredgecolor=C[i], markeredgewidth=1.2, markersize=6,
                label=names[f] + ', leaf-block $L$')
    kk = np.logspace(0, 11, 50)
    ax.plot(kk, kk * U, color=MUTED, lw=0.9, zorder=0)
    ax.annotate(r'$\kappa u$', (kk[-1], kk[-1] * U), xytext=(3, 0), textcoords='offset points', fontsize=7.5, color=INK2, va='center')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(r'$\kappa(HL^{-1})$'); ax.set_ylabel('error vs dense weighted solution')
    ax.legend(fontsize=6.3, loc='upper left')
    fig.tight_layout()
    save(fig, 'm2_weighted')

if 't' in which:
    d1 = load('wt_t1')['t1']
    d2 = load('wt_t2_t3')['t2']
    fig, axs = plt.subplots(1, 3, figsize=(10.2, 3.1))
    names = {'F2': 'F2 wide 1024x2048', 'F2sq': 'F2 square 1536x1536', 'F2tall': 'F2 tall 2048x1024',
             'F4': 'F4 wide conv.', 'F4sq': 'F4 square conv.'}
    for i, f in enumerate(['F2', 'F2sq', 'F2tall', 'F4', 'F4sq']):
        rows = [r for r in d1 if r['family'] == f]
        if not rows:
            continue
        lab = names[f]
        if f in ('F4', 'F4sq'):
            lab = '%s %dx%d' % ('F4 wide' if f == 'F4' else 'F4 square', rows[0]['m'], rows[0]['n'])
        series(axs[0], [r['lam_rel'] for r in rows], [r['err'] for r in rows], i, lab)
    axs[0].set_xscale('log'); axs[0].set_yscale('log'); axs[0].set_xlabel(r'$\lambda/\|H\|_2$')
    axs[0].set_ylabel('error vs SVD filter solution'); axs[0].set_title('(a) accuracy over the $\\lambda$ path', loc='left')
    axs[0].legend(fontsize=6.3, loc='upper right'); axs[0].set_ylim(1e-16, 1e-9)
    wide = [r for r in d2 if r.get('shape', 'wide') == 'wide']
    tall = [r for r in d2 if r.get('shape') == 'tall']
    n = np.array([r['n'] for r in wide])
    series(axs[1], n, [r['t_minnorm'] for r in wide], 0, 'min-norm $H^+b$ (wide, m = n/2)')
    series(axs[1], n, [r['t_tikhonov'] for r in wide], 1, 'Tikhonov (wide, m = n/2)')
    series(axs[1], [r['n'] for r in tall], [r['t_tikhonov'] for r in tall], 2, 'Tikhonov (tall, m = 2n)')
    ref_slope(axs[1], n[3:], wide[3]['t_minnorm'] * 0.5, 1, r'$O(n)$')
    axs[1].set_xscale('log', base=2); axs[1].set_yscale('log'); axs[1].set_xlabel('n (columns)'); axs[1].set_ylabel('seconds (factor + solve)')
    axs[1].set_title('(b) cost', loc='left'); axs[1].legend(fontsize=6.3, loc='upper left')
    series(axs[2], n, [r['err'] for r in wide], 1, 'wide')
    series(axs[2], [r['n'] for r in tall], [r['err'] for r in tall], 2, 'tall')
    axs[2].set_xscale('log', base=2); axs[2].set_yscale('log'); axs[2].set_xlabel('n (columns)')
    axs[2].set_ylabel(r'$\|x-x_\lambda\|/\|x_\lambda\|$'); axs[2].set_ylim(1e-17, 1e-10)
    axs[2].set_title('(c) manufactured $x_\\lambda$, $\\lambda$ = 1e-2', loc='left'); axs[2].legend(fontsize=6.8)
    fig.tight_layout()
    save(fig, 'm2_tikhonov')

if 'bp' in which:
    d = load('apps2_bp')['bp']
    fig, axs = plt.subplots(1, 3, figsize=(10.2, 2.9), gridspec_kw=dict(width_ratios=[1.6, 1, 1]))
    x0 = np.array(d['x0']); n = len(x0); t = np.arange(n)
    ax = axs[0]
    sl = slice(int(0.35 * n), int(0.6 * n))
    ax.plot(t[sl], np.array(d['x_mn'])[sl], color=C[1], lw=0.9, label='minimum-norm $H^+b$')
    ax.plot(t[sl], np.array(d['x_irls'])[sl], color=C[0], lw=1.2, label='basis pursuit (IRLS)')
    ax.plot(t[sl][x0[sl] != 0], x0[sl][x0[sl] != 0], 'o', ms=3.5, color='#52514e', label='true spikes')
    ax.set_xlabel('sample'); ax.set_title('(a) sparse spikes, Ricker wavelet, stride 2', loc='left')
    ax.legend(fontsize=6.5, loc='lower left')
    hi, hd = d['hist']['irls'], d['hist']['dr']
    ax = axs[1]
    series(ax, [h['it'] for h in hi], [h['err'] for h in hi], 0, 'IRLS (weighted min-norm)', markevery=3)
    series(ax, [h['it'] for h in hd], [h['err'] for h in hd], 1, 'Douglas-Rachford', markevery=10)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('iteration'); ax.set_ylabel(r'$\|x_k-x_0\|/\|x_0\|$')
    ax.set_title('(b) convergence', loc='left'); ax.legend(fontsize=6.5)
    ax = axs[2]
    series(ax, [h['t'] for h in hi], [h['err'] for h in hi], 0, 'IRLS: refactor each step', markevery=3)
    series(ax, [h['t'] for h in hd], [h['err'] for h in hd], 1, 'DR: factor once, solve each step', markevery=10)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('seconds'); ax.set_title('(c) wall time', loc='left')
    ax.legend(fontsize=6.5, loc='lower left')
    fig.tight_layout()
    save(fig, 'm2_bp')

if 'db' in which:
    d = load('apps2_db')['db']
    rows = d['rows']
    lam = np.array([r['lam'] for r in rows]); err = np.array([r['err'] for r in rows])
    res = np.array([r['res'] for r in rows]); xn = np.array([r['xnorm'] for r in rows])
    # L-curve corner: maximum curvature of (log res, log xnorm) in the parameter log(lam)
    tt = np.log(lam); r1 = np.gradient(np.log(res), tt); e1 = np.gradient(np.log(xn), tt)
    r2 = np.gradient(r1, tt); e2 = np.gradient(e1, tt)
    curv = (r1 * e2 - r2 * e1) / (r1 ** 2 + e1 ** 2) ** 1.5
    jc = int(np.argmax(curv[1:-1])) + 1
    jd = int(np.argmin(np.abs(lam - d['lam_dp'])))
    fig, axs = plt.subplots(1, 3, figsize=(10.2, 2.9))
    ax = axs[0]
    series(ax, lam, err, 0, None)
    ax.axvline(d['lam_dp'], color=MUTED, lw=0.9)
    ax.axvline(lam[jc], color=MUTED, lw=0.9, ls='--')
    ax.annotate('discrepancy\nprinciple', (d['lam_dp'], err.max()), xytext=(4, -4), textcoords='offset points', fontsize=7, color=INK2, va='top')
    ax.annotate('L-curve\ncorner', (lam[jc], err.max()), xytext=(-4, -4), textcoords='offset points', fontsize=7, color=INK2, va='top', ha='right')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel(r'$\lambda$'); ax.set_ylabel(r'$\|x_\lambda-x_0\|/\|x_0\|$')
    ax.set_title('(a) error along the path', loc='left')
    ax = axs[1]
    series(ax, res, xn, 1, None)
    ax.plot(res[jd], xn[jd], 'o', ms=8, markerfacecolor='none', markeredgecolor=INK2, label='discrepancy, $\\lambda$=%.2g' % lam[jd])
    ax.plot(res[jc], xn[jc], 's', ms=8, markerfacecolor='none', markeredgecolor=INK2, label='max. curvature, $\\lambda$=%.2g' % lam[jc])
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel(r'$\|Hx_\lambda-b\|$'); ax.set_ylabel(r'$\|x_\lambda\|$')
    ax.set_title('(b) L-curve', loc='left'); ax.legend(fontsize=6.5, loc='lower left')
    ax = axs[2]
    t = np.array(d['t'])
    ax.plot(t, d['x0'], color='#b9b8b1', lw=1.4, label='$x_0$')
    ax.plot(t, d['x_dp'], color=C[0], lw=1.0, label='Tikhonov, discrepancy $\\lambda$')
    ax.set_title('(c) reconstruction (1% noise)', loc='left'); ax.legend(fontsize=6.5)
    fig.tight_layout()
    save(fig, 'm2_deblur')
    print('L-curve corner lam %.3g err %.3f; discrepancy lam %.3g err %.3f' % (lam[jc], err[jc], lam[jd], err[jd]))

if 'mri' in which:
    d = load('apps2_mri')['mri']
    t = np.array(d['t'])
    fig, axs = plt.subplots(1, 2, figsize=(7.6, 2.8))
    ax = axs[0]
    ax.plot(t, d['x0'], color='#b9b8b1', lw=1.6, label='phantom')
    ax.plot(t, d['x_mn'], color=C[1], lw=0.9, label='min-norm (err %.2f)' % d['err_minnorm'])
    ax.plot(t, d['x_w'], color=C[0], lw=1.0, label='k-space weighted (err %.2f)' % d['err_weighted'])
    ax.set_xlabel('pixel'); ax.set_title('(a) exact data, m = %d of N = %d' % (d['m'], d['N']), loc='left'); ax.legend(fontsize=6.5)
    ax = axs[1]
    lam = [r['lam'] for r in d['tik']]; e = [r['err'] for r in d['tik']]
    series(ax, lam, e, 0, 'weighted Tikhonov')
    ax.axhline(d['err_minnorm_noisy'], color=C[1], lw=1.0)
    ax.annotate('min-norm, noisy (error %.1e)' % d['err_minnorm_noisy'], (lam[0], d['err_minnorm_noisy']), xytext=(2, -10),
                textcoords='offset points', fontsize=7, color=INK2)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel(r'$\lambda$'); ax.set_ylabel('image error')
    ax.set_title('(b) 2% noise', loc='left'); ax.legend(fontsize=6.5, loc='lower left')
    fig.tight_layout()
    save(fig, 'm2_mri')
print('figures done')
