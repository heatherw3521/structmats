"""Figures, the accuracy table and the numbers quoted in memo 1, from the text files that
hss_examples/minnorm_paper/mp_run_text_results.m writes:

    results/correctness.txt   -> tables/tab_e1.tex
    results/conditioning.txt  -> figures/m1_stability.pdf
    results/scaling.txt       -> figures/m1_complexity.pdf, figures/m1_largescale_accuracy.pdf
    all three                       -> tables/numbers.tex (macros used in the text)

Run from this folder:  python3 make_figures.py   (needs numpy and matplotlib).
MEMO1_RESULTS can point to another results folder; if MEMO1_PNG_DIR is set, PNG previews
of the figures are written there too.
"""
import csv
import math
import os
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.environ.get('MEMO1_RESULTS', os.path.join(HERE, '..', '..', 'hss_examples', 'minnorm_paper', 'results'))
FIG = os.path.join(HERE, 'figures')
TAB = os.path.join(HERE, 'tables')
os.makedirs(FIG, exist_ok=True)
os.makedirs(TAB, exist_ok=True)
U = np.finfo(float).eps / 2

C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
MK = ['o', 's', '^', 'D', 'v', 'P', 'X', '*']
INK, INK2, MUTED, GRID, AXIS = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 9.5, 'axes.labelsize': 9,
    'axes.edgecolor': AXIS, 'axes.labelcolor': INK2, 'axes.linewidth': 0.8, 'axes.grid': True,
    'grid.color': GRID, 'grid.linewidth': 0.6, 'xtick.color': INK2, 'ytick.color': INK2,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 7.5, 'legend.frameon': False,
    'lines.linewidth': 1.6, 'lines.markersize': 5, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.03,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.titlecolor': INK, 'text.color': INK})


def series(ax, x, y, i, label=None, **kw):
    kw.setdefault('marker', MK[i % len(MK)])
    kw.setdefault('markeredgecolor', 'white')
    kw.setdefault('markeredgewidth', 0.8)
    kw.setdefault('color', C[i % len(C)])
    return ax.plot(x, y, label=label, **kw)


def guide(ax, x, y, text, dy=0):
    ax.plot(x, y, color=MUTED, lw=0.9, zorder=0)
    ax.annotate(text, (x[-1], y[-1]), xytext=(3, dy), textcoords='offset points', color=INK2,
                fontsize=7.5, va='center')


def read(name):
    """rows of a results file as dicts (numbers as float), plus its '#' header lines."""
    path = os.path.join(RES, name)
    with open(path) as f:
        lines = f.read().splitlines()
    meta = [l[1:].strip() for l in lines if l.startswith('#')]
    body = [l for l in lines if l and not l.startswith('#')]
    rows = []
    for r in csv.DictReader(body):
        d = {}
        for k, v in r.items():
            try:
                d[k] = float(v)
            except ValueError:
                d[k] = v
        rows.append(d)
    return rows, meta


def e(v, d=1):
    """number in LaTeX scientific notation."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return '--'
    s = '%.*e' % (d, v)
    m, ex = s.split('e')
    return r'$%s{\times}10^{%d}$' % (m, int(ex))


def emath(v, d=1):
    s = '%.*e' % (d, v)
    m, ex = s.split('e')
    return r'%s\times10^{%d}' % (m, int(ex))


PNG = os.environ.get('MEMO1_PNG_DIR')


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + '.pdf'))
    if PNG:
        os.makedirs(PNG, exist_ok=True)
        fig.savefig(os.path.join(PNG, name + '.png'), dpi=160)
    plt.close(fig)


macros = {}
platforms = []

# ---------------------------------------------------------------- accuracy table (tab:e1)
rows, meta = read('correctness.txt')
platforms.append(meta[-1])
L = [r'\begin{tabular}{@{}lrrrrrcccccc@{}}', r'\toprule',
     r'case & $m$ & $n$ & $L$ & rank & $\kappa(H)$ & As.~\ref{as:slack} & residual & solver error & null frac. & vs.\ $A$\\',
     r'\midrule']
def pretty(name):
    name = name.replace('gauss', 'Gaussian').replace('ricker', 'Ricker')
    name = name.replace(' sigma=', r', $\sigma=').replace(' s=', r'$, $s=')
    name = name.replace('m/n = ', r'$m/n=')
    if name.count('$') % 2:
        name += '$'
    return name


def kfmt(v):
    return '%.1f' % v if v < 10 else ('%.0f' % v if v < 100 else e(v))


for r in rows:
    L.append(r'%s & %d & %d & %d & %d & %s & %s & %s & %s & %s & %s\\' % (
        pretty(r['name']), r['m'], r['n'], r['L'], r['maxrank'], kfmt(r['kappa']), 'yes' if r['slack_ok'] == 1 else 'no',
        e(r['res']), e(r['err_H']), e(r['null']), e(r['err_A'])))
L += [r'\bottomrule', r'\end{tabular}']
open(os.path.join(TAB, 'tab_e1.tex'), 'w').write('\n'.join(L) + '\n')
macros['EoneMaxRes'] = emath(max(r['res'] for r in rows))
macros['EoneMaxErr'] = emath(max(r['err_H'] for r in rows))
macros['EoneMaxNull'] = emath(max(r['null'] for r in rows))
big = [r['err_H'] / (r['kappa'] * U) for r in rows]
macros['EoneMaxErrOverKappaU'] = '%d' % math.ceil(max(big))
imax = max(range(len(rows)), key=lambda i: big[i])
small = [r['err_H'] for r in rows if r['kappa'] < 10]
# if the worst ratio comes from a well-conditioned case, say so (its absolute error is tiny)
macros['EoneSmallKappaClause'] = ((r'; the largest ratios occur for $\kappa(H)<10$, where the solver error is at '
                                   r'most $%s$') % emath(max(small))) if rows[imax]['kappa'] < 10 else ''
comp = [r for r in rows if r['family'] in ('F3', 'F5', 'F6')]
macros['EoneCompMin'] = emath(min(r['err_A'] for r in comp))
macros['EoneCompMax'] = emath(max(r['err_A'] for r in comp))
noslack = [r['name'] for r in rows if r['slack_ok'] != 1]
noslack = [pretty(x) for x in noslack]
macros['EoneNoSlack'] = ', '.join(noslack[:-1]) + ' and ' + noslack[-1] if len(noslack) > 1 else ''.join(noslack)
kap = {}
for r in rows:
    kap.setdefault(r['family'], []).append(r['kappa'])

# ---------------------------------------------------------------- Figure: stability
rows, meta = read('conditioning.txt')
platforms.append(meta[-1])
fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.2))
for j, (fam, title) in enumerate([('F2 graded', r'(a) F2$_\alpha$, graded columns, depth 3'),
                                  ('F4 blur', r'(b) F4, Gaussian blur, depth 4')]):
    R = sorted([r for r in rows if r['family'] == fam], key=lambda r: r['kappa'])
    k = np.array([r['kappa'] for r in R])
    ax = axs[j]
    series(ax, k, [max(r['ulv_fwd'], 1e-17) for r in R], 0, 'HSS ULV (this algorithm)')
    series(ax, k, [r['ne_fwd'] for r in R], 1, r'normal equations (Cholesky of $HH^*$)')
    series(ax, k, [r['cg_fwd'] for r in R], 2, r'CGNE (HSS products, $\leq$2000 its)')
    guide(ax, k, k * U, r'$\kappa u$')
    kk = k[k ** 2 * U < 3]
    guide(ax, kk, kk ** 2 * U, r'$\kappa^2 u$')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(1e-17, 5)
    ax.set_xlabel(r'$\kappa(H)$'); ax.set_title(title, loc='left')
    if j == 0:
        ax.set_ylabel('difference from dense QR solution')
    kap[fam] = list(k)
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.03), fontsize=7.5)
fig.tight_layout(rect=(0, 0.07, 1, 1))
save(fig, 'm1_stability')
ulv_ratio = max(r['ulv_fwd'] / (r['kappa'] * U) for r in rows)
macros['EtwoUlvOverKappaU'] = '%d' % math.ceil(ulv_ratio)
macros['EtwoMaxKappaGraded'] = emath(max(kap['F2 graded']))
macros['EtwoMaxKappaBlur'] = emath(max(kap['F4 blur']))
macros['EtwoUlvBwdMax'] = emath(max(r['ulv_bwd'] for r in rows))
macros['EtwoQrBwdMax'] = emath(max(r['qr_bwd'] for r in rows))
ne = [r for r in rows if not math.isnan(r['ne_bwd'])]
macros['EtwoNeBwdMax'] = emath(max(r['ne_bwd'] for r in ne))


def tail(R, bad):
    """smallest kappa from which on every row (R sorted by kappa) is bad."""
    k = None
    for r in reversed(R):
        if not bad(r):
            break
        k = r['kappa']
    return emath(k) if k is not None else '--'


for fam, key in [('F2 graded', 'Graded'), ('F4 blur', 'Blur')]:
    R = sorted([r for r in rows if r['family'] == fam], key=lambda r: r['kappa'])
    macros['EtwoCholFail' + key] = tail(R, lambda r: math.isnan(r['ne_fwd']))
    # CGNE fails: pcg did not converge (flag ~= 0) and the error exceeds 1e-6
    macros['EtwoCgFail' + key] = tail(R, lambda r: r['cg_flag'] != 0 and r['cg_fwd'] > 1e-6)
    good = [r for r in R if r['cg_flag'] == 0]
    macros['EtwoCgItsMin' + key] = '%d' % min(r['cg_its'] for r in good) if good else '--'

# ---------------------------------------------------------------- Figures: complexity, large-n accuracy
rows, meta = read('scaling.txt')
platforms.append(meta[-1])
fams = ['F1', 'F2', 'F3', 'F4', 'F5']
names = {'F1': 'F1', 'F2': 'F2', 'F3': 'F3', 'F4': 'F4', 'F5': 'F5'}


def sel(fam, meth):
    return sorted([r for r in rows if r['family'] == fam and r['method'] == meth], key=lambda r: r['n'])


fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.3))
ax = axs[0]
T = sel('F2', 'ulv'); n = np.array([r['n'] for r in T])
series(ax, n, [r['t'] for r in T], 0, 'HSS ULV, first solve')
series(ax, n, [r['t_repeat'] for r in T], 0, 'HSS ULV, repeat solve', linestyle='--', marker='o', markerfacecolor='white', markeredgecolor=C[0])
D = sel('F2', 'dense_qr')
if D:
    series(ax, [r['n'] for r in D], [r['t'] for r in D], 1, 'dense QR min-norm')
G = sel('F2', 'cgne') + sel('F2', 'cgne_notconverged')
G = sorted(G, key=lambda r: r['n'])
if G:
    series(ax, [r['n'] for r in G], [r['t'] for r in G], 2, 'CGNE (HSS products)')
nn = n[n >= n[len(n) // 2]]
guide(ax, nn, T[-1]['t'] * 1.6 * nn / nn[-1], r'$O(n)$', dy=2)
if D:
    nd = np.array([r['n'] for r in D], dtype=float)
    guide(ax, nd, D[-1]['t'] * 1.6 * (nd / nd[-1]) ** 3, r'$O(n^3)$')
ax.set_xscale('log', base=2); ax.set_yscale('log'); ax.set_xlabel('n (columns)'); ax.set_ylabel('seconds')
ax.set_title(r'(a) F2, $m=n/2$', loc='left'); ax.legend(loc='lower right', fontsize=6.8)
ax = axs[1]
for i, f in enumerate(fams):
    T = sel(f, 'ulv')
    if not T:
        continue
    series(ax, [r['n'] for r in T], [1e6 * r['t'] / r['n'] for r in T], i, names[f])
ax.set_xscale('log', base=2); ax.set_yscale('log'); ax.set_xlabel('n (columns)')
ax.set_ylabel(r'first-solve time / $n$  ($\mu$s)'); ax.set_title('(b) time per column', loc='left')
per = {f: np.log10([1e6 * r['t'] / r['n'] for r in sel(f, 'ulv')]) for f in fams if sel(f, 'ulv')}
lo = min(v.min() for v in per.values()); hi = max(v.max() for v in per.values())
ylo, yhi = lo - 0.12, hi + 0.12
ax.set_ylim(10 ** ylo, 10 ** yhi)
ax.set_yticks([t for t in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000) if 10 ** ylo <= t <= 10 ** yhi])
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%g'))
ax.yaxis.set_minor_formatter(mticker.NullFormatter())
# legend in the widest empty band between the curves (fallback: best)
band = sorted((v.min(), v.max()) for v in per.values())
gaps = [(band[k + 1][0] - max(b[1] for b in band[:k + 1]), k) for k in range(len(band) - 1)]
g, k = max(gaps) if gaps else (0, 0)
if g > 0.25:
    yc = (max(b[1] for b in band[:k + 1]) + g / 2 - ylo) / (yhi - ylo)
    ax.legend(loc='center', bbox_to_anchor=(0.5, yc), ncol=5, fontsize=6.8, handlelength=1.6, columnspacing=1.0)
else:
    ax.legend(loc='best', fontsize=6.8, ncol=2)
fig.tight_layout()
save(fig, 'm1_complexity')


def slope(T, key='t', last=4):
    T = [r for r in T if r[key] > 0][-last:]
    x = np.log([r['n'] for r in T]); y = np.log([r[key] for r in T])
    return np.polyfit(x, y, 1)[0]


sl = [slope(sel(f, 'ulv')) for f in fams if sel(f, 'ulv')]
macros['EfourSlopeMin'] = '%.2f' % min(sl)
macros['EfourSlopeMax'] = '%.2f' % max(sl)
slr = [slope(sel(f, 'ulv'), 't_repeat') for f in fams if sel(f, 'ulv')]
macros['EfourSlopeRepeatMin'] = '%.2f' % min(slr)
macros['EfourSlopeRepeatMax'] = '%.2f' % max(slr)
last = [sel(f, 'ulv')[-1] for f in fams if sel(f, 'ulv')]
rr = sorted(round(r['t'] / r['t_repeat']) for r in last)          # first / repeat solve, largest n
macros['EfourRatioRepeat'] = '%d' % rr[0] if rr[0] == rr[-1] else '%d--%d' % (rr[0], rr[-1])
rm = sorted(round(r['t_repeat'] / r['t_matvec']) for r in last)   # repeat solve in HSS products
macros['EfourRepeatMatvec'] = '%d' % rm[0] if rm[0] == rm[-1] else '%d--%d' % (rm[0], rm[-1])
T = sel('F2', 'ulv')
macros['EfourMaxN'] = '2^{%d}' % round(math.log2(max(r['n'] for r in rows)))
if D:
    dl = D[-1]; tl = [r for r in T if r['n'] == dl['n']][0]
    macros['EfourDenseN'] = '%d' % dl['n']
    macros['EfourDenseRatio'] = '%.0f' % (dl['t'] / tl['t'])
    macros['EfourDenseSlope'] = '%.1f' % slope(D)
if G:
    gl = G[-1]; tl = [r for r in T if r['n'] == gl['n']][0]
    macros['EfourCgN'] = '%d' % gl['n']
    macros['EfourCgRatio'] = '%.0f' % (gl['t'] / tl['t'])
    v = gl['t'] / tl['t_repeat']                                     # CGNE vs repeat solve, 2 digits
    macros['EfourCgRatioRepeat'] = '%d' % round(v, 1 - int(math.floor(math.log10(v))))
    macros['EfourCgIts'] = '%d--%d' % (min(r['iters'] for r in G), max(r['iters'] for r in G))

fig, axs = plt.subplots(1, 2, figsize=(8.0, 3.0))
for i, f in enumerate(fams):
    T = sel(f, 'ulv')
    if not T:
        continue
    series(axs[0], [r['n'] for r in T], [r['err'] for r in T], i, names[f])
    series(axs[1], [r['n'] for r in T], [r['res'] for r in T], i)
axs[0].set_ylabel(r'$\|x-x_\star\|/\|x_\star\|$'); axs[1].set_ylabel(r'$\|Hx-b\|/\|b\|$')
axs[0].set_title('(a) error vs. manufactured solution', loc='left'); axs[1].set_title('(b) residual', loc='left')
for ax in axs:
    ax.set_xscale('log', base=2); ax.set_yscale('log'); ax.set_xlabel('n (columns)'); ax.set_ylim(1e-17, 1e-8)
axs[0].legend(loc='upper left', fontsize=6.8, ncol=2)
fig.tight_layout()
save(fig, 'm1_largescale_accuracy')
ok = [r for r in rows if r['method'] == 'ulv']
macros['EthreeMaxRes'] = emath(max(r['res'] for r in ok))
nf1 = [r for r in ok if r['family'] != 'F1']
macros['EthreeErrMin'] = emath(min(r['err'] for r in nf1))
macros['EthreeErrMax'] = emath(max(r['err'] for r in nf1))
F1 = sel('F1', 'ulv')
macros['EthreeFoneFirst'] = emath(F1[0]['err']); macros['EthreeFoneLast'] = emath(F1[-1]['err'])
macros['EthreeFoneN'] = '2^{%d}' % round(math.log2(F1[-1]['n']))

# ---------------------------------------------------------------- kappa column of the families table
def fmt0(v):
    """1.3e3 -> 1{\times}10^{3} style, small numbers plain."""
    if v < 10:
        return '%.1f' % v
    if v < 100:
        return '%.0f' % v
    ex = int(math.floor(math.log10(v)))
    m = round(v / 10 ** ex)
    if m == 10:
        m, ex = 1, ex + 1
    return r'10^{%d}' % ex if m == 1 else r'%d{\times}10^{%d}' % (m, ex)


def krange(v):
    lo, hi = min(v), max(v)
    if hi / lo < 3:
        return r'$\approx%s$' % fmt0(lo)
    return r'$%s$--$%s$' % (fmt0(lo), fmt0(hi))


for fam, key in [('F1', 'KFone'), ('F2', 'KFtwo'), ('F3', 'KFthree'), ('F5', 'KFfive'), ('F6', 'KFsix')]:
    macros[key] = krange(kap[fam])
macros['KFtwoG'] = krange(kap['F2 graded'])
macros['KFfour'] = krange(kap['F4'] + kap['F4 blur'])


def platform(line):
    """'<version> | <computer> | threads <k> | <date>' -> (name, machine, threads)."""
    parts = [p.strip() for p in line.split('|')]
    m = re.search(r'\((R\d{4}[ab])\)', parts[0])
    name = 'MATLAB ' + m.group(1) if m else 'GNU Octave ' + parts[0]
    machine = {'MACA64': 'Apple silicon Mac', 'MACI64': 'Intel Mac', 'GLNXA64': 'Linux',
               'PCWIN64': 'Windows'}.get(parts[1], 'Linux' if 'linux' in parts[1] else parts[1])
    threads = parts[2].replace('threads', '').strip() if len(parts) > 2 else '?'
    return name, machine, threads


acc = sorted(set(platform(p)[0] for p in platforms[:2]))     # correctness, conditioning
tim = platform(platforms[2])                                  # scaling (timings, large n)
macros['TimingPlatform'] = '%s (%s, %s threads)' % tim
if acc == [tim[0]]:
    macros['PlatformSentence'] = r'All results were computed with \TimingPlatform.'
else:
    macros['PlatformSentence'] = (r'The timings and the large-$n$ results were computed with \TimingPlatform; the '
                                  r'accuracy results of Sections~\ref{subsec:correct} and~\ref{subsec:stab} with '
                                  + ' and '.join(acc) + ', from the same files.')

with open(os.path.join(TAB, 'numbers.tex'), 'w') as f:
    f.write('%% generated by make_figures.py from %s\n' % os.path.relpath(RES, HERE))
    for p in platforms:
        f.write('%% %s\n' % p)
    for k, v in macros.items():
        f.write('\\newcommand{\\%s}{%s}\n' % (k, v))
for k, v in macros.items():
    print('%-24s %s' % (k, v))
