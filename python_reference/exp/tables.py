"""LaTeX table snippets generated from results (no hand transcription)."""
import os, json
import numpy as np
from common import load
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get('MEMO1_TABLES', os.path.join(HERE, '..', '..', 'notes', 'minnorm_v2', 'tables'))
# Octave cross-check output (one row per matrix: name m n L slack_ok nbad residual
# e_vs_collapse e_vs_relax e_prototype e_weighted e_tikhonov), current and legacy solver
XVAL_CURRENT = os.environ.get('XVAL_CURRENT', os.path.join(HERE, 'results', 'xval2_results_current.txt'))
XVAL_LEGACY = os.environ.get('XVAL_LEGACY', os.path.join(HERE, 'results', 'xval2_results_legacy.txt'))
os.makedirs(OUT, exist_ok=True)

def e(v, d=1):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return '--'
    s = '%.*e' % (d, v)
    m, ex = s.split('e')
    return r'$%s{\times}10^{%d}$' % (m, int(ex))

def e1():
    rows = load('m1_correctness')['rows']
    L = [r'\begin{tabular}{@{}lrrrrrrrrc@{}}', r'\toprule',
         r'case & $m$ & $n$ & $L$ & rank & $\kappa(H)$ & residual & solver error & null frac. & vs.\ $A$\\', r'\midrule']
    for r in rows:
        x = r['relax']
        nm = r['name'].replace('->', r'$\to$').replace('_', r'\_')
        L.append('%s & %d & %d & %d & %d & %s & %s & %s & %s & %s\\\\' % (
            nm, r['m'], r['n'], r['L'], r['maxrank'], e(r['kappa']), e(x['res']), e(x['err_qr']), e(x['null']),
            e(r.get('err_vs_A')) if 'err_vs_A' in r else '--'))
    L += [r'\bottomrule', r'\end{tabular}']
    open(os.path.join(OUT, 'tab_e1.tex'), 'w').write('\n'.join(L) + '\n')

def scaling():
    from common import fit_slope
    d = load('m3_scaling')
    names = {'F1': 'F1 LR+block-diag.', 'F2': 'F2 random HSS', 'F3': 'F3 Cauchy (proxy)', 'F4': 'F4 conv. (proxy)', 'F5': 'F5 NUDFT (proxy)'}
    L = [r'\begin{tabular}{@{}lrrrrrrr@{}}', r'\toprule',
         r'family & largest $n$ & rank & build (s) & factor (s) & solve (s) & slope & max error\\', r'\midrule']
    for f in ['F1', 'F2', 'F3', 'F4', 'F5']:
        R = d[f]; n = [r['n'] for r in R]; tf = [r['t_factor'] for r in R]
        b = R[-1].get('build')
        L.append('%s & %d & %d & %s & %.2f & %.2f & %.2f & %s\\\\' % (names[f], n[-1], R[-1]['maxrank'],
                 ('%.1f' % b) if b else '--', tf[-1], R[-1]['t_solve'], fit_slope(n, tf, 4), e(max(r['err'] for r in R))))
    L += [r'\bottomrule', r'\end{tabular}']
    open(os.path.join(OUT, 'tab_scaling.tex'), 'w').write('\n'.join(L) + '\n')

def xval():
    def read(p):
        return {r[0]: r for r in (l.split() for l in open(p) if l.strip())}
    cur, leg = read(XVAL_CURRENT), read(XVAL_LEGACY)
    def f(v):
        v = float(v)
        return '--' if np.isnan(v) else e(v)
    L = [r'\begin{tabular}{@{}lrrcccccc@{}}', r'\toprule',
         r'matrix & $m\times n$ & $L$ & As.~\ref{as:slack} & residual & current vs.\ Python & legacy vs.\ Python & weighted & Tikhonov\\', r'\midrule']
    for name, r in cur.items():
        _, m, n, Lv, ok, nbad, res, emn, emnr, erel, ew, et = r
        # current: H\b (square: minnormfactor) against the port's relaxed rule;
        # legacy: legacy H\b against the port's slack='collapse' rule
        L.append(r'%s & %s$\times$%s & %s & %s & %s & %s & %s & %s & %s\\' % (
            name.replace('_', r'\_'), m, n, Lv, 'yes' if ok == '1' else 'no', f(res), f(emnr),
            f(leg[name][7]), f(ew), f(et)))
    L += [r'\bottomrule', r'\end{tabular}']
    open(os.path.join(OUT, 'tab_xval.tex'), 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    e1(); scaling(); xval()
    for t in ['tab_e1', 'tab_scaling', 'tab_xval']:
        print(open(os.path.join(OUT, t + '.tex')).read())
