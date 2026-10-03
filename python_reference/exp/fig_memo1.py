"""Figures for memo 1 from results/*.json."""
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from common import style, series, ref_slope, save, load, C, MK, INK2, MUTED, fit_slope

style()
U = np.finfo(float).eps / 2
which = sys.argv[1:] or ['stab', 'scale', 'sweeps', 'a1', 'a2', 'a34']

# ---------------------------------------------------------------- E2
if 'stab' in which:
    d = load('m2_conditioning')
    d2 = load('m2b_deep')
    fig, axs = plt.subplots(2, 3, figsize=(10.2, 5.4), sharex='col')
    meth = [('ulv', 0, 'HSS ULV, factored root solve (current code)'), ('qr', 1, 'dense QR of $H^*$'),
            ('normal_eq', 2, 'normal eq. (Cholesky $HH^*$)'), ('cgne', 3, 'CGNE (2000 its max)'),
            ('ulv_pinv', 7, 'HSS ULV, explicit pinv at root (legacy code)')]
    for j, (fam, title) in enumerate([('graded', '(a) F2 with column grading, depth 3'),
                                      ('blur', '(b) Gaussian deconvolution (F4), depth 4')]):
        rows = d[fam]
        k = np.array([r['kappa'] for r in rows])
        o = np.argsort(k); k = k[o]
        for row, key, lab in [(0, 'fwd', 'forward error (60-digit ref.)'), (1, 'bwd', 'backward error')]:
            ax = axs[row, j]
            for mk, ci, ml in meth:
                if row == 0 and mk == 'ulv_pinv':
                    continue          # same forward error as 'ulv'
                v = np.array([r[mk][key] for r in rows], dtype=float)[o]
                series(ax, k, v, ci, ml)
            if key == 'fwd':
                ax.plot(k, k * U, color=MUTED, lw=0.9, zorder=0)
                ax.plot(k, np.minimum(k ** 2 * U, 10), color=MUTED, lw=0.9, zorder=0)
                ax.annotate(r'$\kappa u$', (k[-1], k[-1] * U), xytext=(3, 0), textcoords='offset points', fontsize=7.5, color=INK2, va='center')
                kk = k[k ** 2 * U < 1]
                ax.annotate(r'$\kappa^2 u$', (kk[-1], kk[-1] ** 2 * U), xytext=(3, -2), textcoords='offset points', fontsize=7.5, color=INK2, va='center')
                ax.set_ylim(1e-17, 3)
            else:
                ax.plot(k, np.full_like(k, U), color=MUTED, lw=0.9, zorder=0)
                ax.annotate(r'$u$', (k[-1], U), xytext=(4, 0), textcoords='offset points', fontsize=7.5, color=INK2, va='center')
                ax.set_ylim(1e-18, 3)
            ax.set_xscale('log'); ax.set_yscale('log')
            if j == 0:
                ax.set_ylabel(lab)
            if row == 0:
                ax.set_title(title, loc='left')
            if row == 1:
                ax.set_xlabel(r'$\kappa(H)$')
    # (c) deeper tree, dense-QR reference
    rows = d2['rows']
    k = np.array([r['kappa'] for r in rows])
    ax = axs[0, 2]
    series(ax, k, [r['ulv']['fwd'] for r in rows], 0, None)
    ax.plot(k, k * U, color=MUTED, lw=0.9, zorder=0)
    ax.annotate(r'$\kappa u$', (k[-1], k[-1] * U), xytext=(3, 0), textcoords='offset points', fontsize=7.5, color=INK2, va='center')
    ax.set_ylim(1e-17, 3); ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_ylabel('difference from dense QR')
    ax.set_title('(c) F2 graded, %dx%d, depth %d' % (d2['m'], d2['n'], d2['L']), loc='left')
    ax = axs[1, 2]
    series(ax, k, [r['ulv']['bwd'] for r in rows], 0, None)
    series(ax, k, [r['qr_bwd'] for r in rows], 1, None)
    series(ax, k, [r['ulv_pinv']['bwd'] for r in rows], 7, None)
    ax.plot(k, np.full_like(k, U), color=MUTED, lw=0.9, zorder=0)
    ax.annotate(r'$u$', (k[-1], U), xytext=(4, 0), textcoords='offset points', fontsize=7.5, color=INK2, va='center')
    ax.set_ylim(1e-18, 3); ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel(r'$\kappa(H)$')
    h, l = axs[1, 0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=5, bbox_to_anchor=(0.5, -0.01), fontsize=7.3)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, 'm1_stability')

# ---------------------------------------------------------------- E3/E4
if 'scale' in which:
    d = load('m3_scaling')
    names = {'F1': 'F1 LR+block-diag. (k=3, bs 16)', 'F2': 'F2 random HSS (k=10)', 'F3': 'F3 interlaced Cauchy',
             'F4': 'F4 decimated Gaussian conv.', 'F5': 'F5 NUDFT Cauchy-like'}
    fams = [f for f in ['F1', 'F2', 'F3', 'F4', 'F5'] if f in d]
    fig, axs = plt.subplots(1, 3, figsize=(10.2, 3.3))
    for i, f in enumerate(fams):
        n = np.array([r['n'] for r in d[f]]); tf = np.array([r['t_factor'] for r in d[f]])
        ts = np.array([r['t_solve'] for r in d[f]])
        series(axs[0], n, tf, i, names[f])
        series(axs[1], n, ts, i, names[f])
    n2 = np.array([r['n'] for r in d['F2']])
    ref_slope(axs[0], n2[3:], 0.5 * d['F2'][3]['t_factor'], 1, r'$O(n)$')
    ref_slope(axs[1], n2[3:], 0.5 * d['F2'][3]['t_solve'], 1, r'$O(n)$')
    for ax, t in [(axs[0], '(a) factorization (all levels)'), (axs[1], '(b) solve with stored factors')]:
        ax.set_xscale('log', base=2); ax.set_yscale('log'); ax.set_xlabel('n (columns)'); ax.set_ylabel('seconds')
        ax.set_title(t, loc='left')
    axs[0].legend(loc='upper left', fontsize=6.8)
    if 'dense_F2' in d:
        D = d['dense_F2']
        n = np.array([r['n'] for r in D])
        series(axs[2], n, [r['t_dense'] for r in D], 0, 'dense min-norm (QR)')
        series(axs[2], n, [r['t_factor'] + r['t_solve'] for r in D], 1, 'HSS factor + solve')
        series(axs[2], n, [r['t_cgne'] for r in D], 2, 'CGNE to 1e-12 (HSS products)')
        axs[2].set_xscale('log', base=2); axs[2].set_yscale('log'); axs[2].set_xlabel('n (F2, m = n/2)')
        axs[2].set_ylabel('seconds'); axs[2].set_title('(c) against dense and iterative', loc='left')
        axs[2].legend(loc='upper left', fontsize=6.8)
    fig.tight_layout()
    save(fig, 'm1_timing')
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0))
    for i, f in enumerate(fams):
        n = np.array([r['n'] for r in d[f]])
        series(axs[0], n, [r['err'] for r in d[f]], i, names[f])
        series(axs[1], n, [r['res'] for r in d[f]], i)
    axs[0].set_ylabel(r'$\|x-x_\star\|/\|x_\star\|$'); axs[1].set_ylabel(r'$\|Hx-b\|/\|b\|$')
    axs[0].set_title('(a) error vs manufactured min-norm solution', loc='left')
    axs[1].set_title('(b) residual', loc='left')
    for ax in axs:
        ax.set_xscale('log', base=2); ax.set_yscale('log'); ax.set_xlabel('n (columns)'); ax.set_ylim(1e-17, 1e-8)
    axs[0].legend(loc='upper left', fontsize=6.8)
    fig.tight_layout()
    save(fig, 'm1_largescale_accuracy')

# ---------------------------------------------------------------- E5/E6
if 'sweeps' in which:
    d = load('m5_sweeps')
    fig, axs = plt.subplots(1, 3, figsize=(10.2, 3.2))
    for ci, bsz in [(0, 64), (2, 128)]:
        rows = [r for r in d['rank'] if r['bs'] == bsz]
        k = np.array([r['k'] for r in rows]); t = np.array([r['t_factor'] for r in rows])
        fl = np.array([r['flops'] for r in rows])
        sc = np.exp(np.mean(np.log(t / fl)))
        n0 = rows[0]['n0']
        series(axs[0], k, t, ci, 'leaves %dx%d' % (n0, 2 * n0))
        axs[0].plot(k, sc * fl, color=C[ci], lw=0.9, ls='--', zorder=0)
    axs[0].plot([], [], color=MUTED, lw=0.9, ls='--', label='flop count (Sec. 5), scaled')
    axs[0].set_xscale('log', base=2); axs[0].set_yscale('log'); axs[0].set_xlabel('HSS rank k'); axs[0].set_ylabel('factor time (s)')
    axs[0].set_title('(a) rank, n = 65 536, fixed leaf size', loc='left')
    axs[0].legend(loc='upper left', fontsize=6.8)
    bs = np.array([r['bs'] for r in d['blocksize']]); t = np.array([r['t_factor'] for r in d['blocksize']])
    series(axs[1], bs, t, 1)
    axs[1].set_xscale('log', base=2); axs[1].set_yscale('log'); axs[1].set_xlabel('blocksize'); axs[1].set_ylabel('factor time (s)')
    axs[1].set_title('(b) blocksize, n = 65 536, k = 10', loc='left')
    A = d['aspect']
    Ac = [r for r in A if r.get('collapse')]
    rr = np.array([r['ratio'] for r in A]); rc = np.array([r['ratio'] for r in Ac])
    tr = np.array([r['relax']['t_factor'] for r in A]); tc = np.array([r['collapse']['t_factor'] for r in Ac])
    series(axs[2], rc, tc, 7, "legacy rule (simulated): merge leaf level")
    series(axs[2], rr, tr, 0, r"relaxed $p'=\min(p,\ell+n)$ (current code)")
    for r_, t_, a in zip(rc, tc, Ac):
        c = a['collapse']['collapses']
        if c:
            axs[2].annotate('%d' % c, (r_, t_), xytext=(0, 5), textcoords='offset points', fontsize=7, color=INK2, ha='center')
    axs[2].set_yscale('log'); axs[2].set_xlabel('aspect ratio m/n'); axs[2].set_ylabel('factor time (s)')
    axs[2].set_title('(c) aspect ratio, n = 16 384, k = 16', loc='left')
    axs[2].legend(loc='upper left', fontsize=6.8)
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    for ax_, tk in [(axs[0], [0.3, 0.5, 1, 2, 3]), (axs[1], [0.5, 1, 2, 5]), (axs[2], [0.1, 0.2, 0.5, 1])]:
        ax_.yaxis.set_major_locator(FixedLocator(tk)); ax_.yaxis.set_minor_locator(NullLocator())
        ax_.yaxis.set_major_formatter(FuncFormatter(lambda v, _: ('%g' % v)))
    fig.tight_layout()
    save(fig, 'm1_sweeps')

# ---------------------------------------------------------------- A1
if 'a1' in which:
    d = load('a1_nudft')
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 2.8), gridspec_kw=dict(width_ratios=[2.2, 1]))
    g = d['gap']
    tt = np.array(g['tt']); ax = axs[0]
    ax.plot(tt, g['ftrue'], color='#b9b8b1', lw=1.4, label='true band-limited signal')
    ax.plot(tt, g['fmn'], color=C[0], lw=1.1, label='minimum-norm interpolant')
    ax.plot(g['xs'], g['fx'], 'o', ms=2.2, color=C[1], label='samples (m = %d, N = %d)' % (g['m'], g['N']))
    ax.axvspan(0.42, 0.52, color='#f0efec', zorder=0)
    ax.set_xlim(0.3, 0.64); ax.set_xlabel('x'); ax.set_title('(a) gap in the sampling (shaded)', loc='left')
    yl = ax.get_ylim(); ax.set_ylim(yl[0], yl[1] + 0.45 * (yl[1] - yl[0]))
    ax.legend(loc='upper left', fontsize=6.8, ncol=3)
    ax = axs[1]
    labels = ['jitter', 'gap']
    xi = np.arange(2)
    ax.bar(xi - 0.17, [d[c]['trees']['index_tree']['maxrank'] for c in labels], 0.3, color=C[0], label='index tree (hss_constructor)')
    ax.bar(xi + 0.17, [d[c]['trees']['geometry_tree']['maxrank'] for c in labels], 0.3, color=C[1], label='geometry-aligned tree')
    for j, c in enumerate(labels):
        for off, key in [(-0.17, 'index_tree'), (0.17, 'geometry_tree')]:
            v = d[c]['trees'][key]['maxrank']
            ax.annotate(str(v), (j + off, v), xytext=(0, 2), textcoords='offset points', ha='center', fontsize=7, color=INK2)
    ax.set_xticks(xi); ax.set_xticklabels(['jittered', 'with gap']); ax.set_ylabel('max HSS rank')
    ax.set_title('(b) tree vs sampling', loc='left'); ax.legend(loc='upper left', fontsize=6.5)
    ax.set_ylim(0, 85)
    fig.tight_layout()
    save(fig, 'm1_a1_nudft')

# ---------------------------------------------------------------- A2
if 'a2' in which:
    d = load('a2_a3_a4_a2')['a2']
    t = np.array(d['t'])
    fig, axs = plt.subplots(1, 3, figsize=(10.2, 2.7))
    ax = axs[0]
    ax.plot(t, d['x0'], color='#b9b8b1', lw=1.4, label='true signal $x_0$')
    ax.plot(t, d['x_mn'], color=C[0], lw=1.0, label='min-norm, exact data')
    ax.set_title('(a) exact data: $x^\\star=P_{\\mathrm{Range}(H^*)}x_0$', loc='left'); ax.legend(fontsize=6.8)
    ax = axs[1]
    ax.plot(t, d['x_noisy'], color=C[7], lw=0.6, label='min-norm, 0.1% noise')
    ax.plot(t, d['x0'], color='#b9b8b1', lw=1.2, label='$x_0$')
    ax.set_title('(b) noisy data (motivates memo 2)', loc='left'); ax.legend(fontsize=6.8)
    ax = axs[2]
    sv = np.array(d['sv'])
    ax.semilogy(np.arange(1, len(sv) + 1), sv, color=C[0], lw=1.4)
    ax.set_xlabel('index'); ax.set_ylabel('singular value'); ax.set_title('(c) spectrum, $\\kappa$ = %.1e' % d['kappa'], loc='left')
    fig.tight_layout()
    save(fig, 'm1_a2_deconv')

# ---------------------------------------------------------------- A3 / A4
if 'a34' in which:
    d3 = load('a2_a3_a4_a3')['a3']
    d4 = load('a2_a3_a4_a4')['a4']
    fig, axs = plt.subplots(1, 3, figsize=(10.2, 3.0))
    n = np.array([r['n'] for r in d3])
    series(axs[0], n, [r['rank_C'] for r in d3], 0, 'Cauchy-like $C=F_mTD^*F_n^*$')
    nT = [r['n'] for r in d3 if 'rank_T' in r]
    series(axs[0], nT, [r['rank_T'] for r in d3 if 'rank_T' in r], 1, 'Toeplitz $T$ itself')
    axs[0].set_xscale('log', base=2); axs[0].set_yscale('log'); axs[0].set_xlabel('n (m = n/2)'); axs[0].set_ylabel('max HSS rank (tol 1e-12)')
    axs[0].set_title('(a) A3: random Toeplitz', loc='left'); axs[0].legend(fontsize=6.8)
    series(axs[1], n, [r['t_factor'] for r in d3], 0, 'factor')
    series(axs[1], n, [r['t_solve_incl_fft'] for r in d3], 1, 'solve incl. FFTs')
    ref_slope(axs[1], n[2:], d3[2]['t_factor'] * 0.5, 1, r'$O(n)$')
    axs[1].set_xscale('log', base=2); axs[1].set_yscale('log'); axs[1].set_xlabel('n'); axs[1].set_ylabel('seconds')
    axs[1].set_title('(b) A3: cost', loc='left'); axs[1].legend(fontsize=6.8)
    for i, order in enumerate(['lex', 'morton']):
        rows = [r for r in d4 if r['order'] == order]
        N = np.array([r['N'] for r in rows])
        series(axs[2], N, [r['maxrank'] for r in rows], i, {'lex': 'lexicographic order', 'morton': 'Morton (Z) order'}[order])
    N = np.array([r['N'] for r in d4 if r['order'] == 'morton'])
    r0 = [r['maxrank'] for r in d4 if r['order'] == 'morton'][0]
    ref_slope(axs[2], N, r0 * 0.7, 0.5, r'$O(\sqrt{N})$')
    axs[2].set_xscale('log', base=2); axs[2].set_yscale('log'); axs[2].set_xlabel('N = pixels')
    axs[2].set_ylabel('max HSS rank'); axs[2].set_title('(c) A4: 2-D blur', loc='left'); axs[2].legend(fontsize=6.8)
    fig.tight_layout()
    save(fig, 'm1_a3_a4')
print('figures done')
