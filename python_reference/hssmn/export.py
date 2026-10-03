import numpy as np
import scipy.io as sio


def _cell(lst):
    c = np.empty((len(lst), 1), dtype=object)
    for i, v in enumerate(lst):
        c[i, 0] = v
    return c


def generators_struct(H, blocksize=None):
    L = H.L
    G = {'L': float(L)}
    G['rb'] = _cell([np.asarray(H.rb[l], dtype=float)[None, :] for l in range(L + 1)])
    G['cb'] = _cell([np.asarray(H.cb[l], dtype=float)[None, :] for l in range(L + 1)])
    G['D'] = _cell(list(H.D))
    U = [np.zeros((0, 0))] + [_cell(list(H.U[l])) for l in range(1, L + 1)]
    V = [np.zeros((0, 0))] + [_cell(list(H.V[l])) for l in range(1, L + 1)]
    G['U'] = _cell(U)
    G['V'] = _cell(V)
    G['B12'] = _cell([_cell([b[0] for b in H.B[l]]) for l in range(L)])
    G['B21'] = _cell([_cell([b[1] for b in H.B[l]]) for l in range(L)])
    if blocksize is not None:
        G['blocksize'] = float(blocksize)
    return G


def save_generators(fname, H, extra=None, blocksize=None):
    d = {'G': generators_struct(H, blocksize)}
    if extra:
        d.update(extra)
    sio.savemat(fname, d, do_compression=True)
