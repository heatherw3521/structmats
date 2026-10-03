"""Shared helpers for the experiment scripts: timing, results I/O, figure style."""
import json
import os
import platform
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(ROOT, 'results')
FIG = os.path.join(ROOT, 'figures')
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

# validated categorical palette (dataviz reference instance, light mode, fixed order)
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
MK = ['o', 's', '^', 'D', 'v', 'P', 'X', '*']
INK, INK2, MUTED, GRID, AXIS = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'


def style():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans',
        'font.size': 9,
        'axes.titlesize': 9.5,
        'axes.labelsize': 9,
        'axes.edgecolor': AXIS,
        'axes.labelcolor': INK2,
        'axes.linewidth': 0.8,
        'axes.grid': True,
        'grid.color': GRID,
        'grid.linewidth': 0.6,
        'grid.linestyle': '-',
        'xtick.color': INK2,
        'ytick.color': INK2,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 7.5,
        'legend.frameon': False,
        'lines.linewidth': 1.6,
        'lines.markersize': 5,
        'lines.solid_capstyle': 'round',
        'lines.solid_joinstyle': 'round',
        'figure.dpi': 150,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.03,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.titlecolor': INK,
        'text.color': INK,
    })


def series(ax, x, y, i, label=None, **kw):
    """line + marker with a white ring (secondary encoding = marker shape)."""
    kw.setdefault('marker', MK[i % len(MK)])
    kw.setdefault('markeredgecolor', 'white')
    kw.setdefault('markeredgewidth', 0.8)
    kw.setdefault('color', C[i % len(C)])
    return ax.plot(x, y, label=label, **kw)


def ref_slope(ax, x, y0, p, text, xpos=None):
    """thin muted reference line y = y0 (x/x[0])^p with a direct label."""
    x = np.asarray(x, dtype=float)
    y = y0 * (x / x[0]) ** p
    ax.plot(x, y, color=MUTED, lw=0.9, zorder=0)
    xi = x[-1] if xpos is None else xpos
    yi = y0 * (xi / x[0]) ** p
    ax.annotate(text, (xi, yi), xytext=(3, 0), textcoords='offset points', color=INK2,
                fontsize=7.5, va='center')


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + '.pdf'))
    fig.savefig(os.path.join(FIG, name + '.png'), dpi=160)
    plt.close(fig)


def timeit(fn, reps=5, warmup=1):
    import gc
    out = None
    for _ in range(warmup):
        out = None
        gc.collect()
        out = fn()
    ts = []
    for _ in range(reps):
        out = None          # free the previous result before the next run (memory)
        gc.collect()
        t0 = time.perf_counter()
        out = fn()
        ts.append(time.perf_counter() - t0)
    return float(np.median(ts)), ts, out


def machine():
    import scipy
    return dict(platform=platform.platform(), python=platform.python_version(),
                numpy=np.__version__, scipy=scipy.__version__,
                cpus=os.cpu_count(), blas_threads=os.environ.get('OPENBLAS_NUM_THREADS', 'default'))


def dump(name, obj):
    obj = dict(obj)
    obj['_machine'] = machine()
    obj['_time'] = time.strftime('%Y-%m-%d %H:%M:%S')
    with open(os.path.join(RES, name + '.json'), 'w') as f:
        json.dump(obj, f, indent=1, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o))


def load(name):
    with open(os.path.join(RES, name + '.json')) as f:
        return json.load(f)


def fit_slope(N, t, tail=3):
    N = np.asarray(N, float)[-tail:]
    t = np.asarray(t, float)[-tail:]
    return float(np.polyfit(np.log(N), np.log(t), 1)[0])
