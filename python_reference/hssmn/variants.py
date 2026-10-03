"""
Weighted minimum-norm and Tikhonov solves built on the unweighted HSS
minimum-norm solver (memo 2).

  weighted:  min ||L x||_2  s.t.  H x = b,   L = diag(w) or blkdiag(L_tau)
             (blocks conforming to the leaf column partition)
             -> y = L x,  (H L^{-1}) y = b,  x = L^{-1} y
  Tikhonov:  min ||S (H x - b)||^2 + lam^2 ||L x||^2   (S, L diagonal / leaf-block)
             -> standard form  Ht = S H L^{-1}, bt = S b
             -> augmented min-norm problem  [Ht, lam I] [y; s] = bt
                (interleaved per leaf so that it is HSS on the same tree)
"""
import numpy as np
from .hss import HSS
from .ulv import MinNormFactor


def right_blockdiag(H, Minv_blocks):
    """H blkdiag(M_tau) for leaf-conforming square blocks M_tau (given as the
    matrices to right-multiply by)."""
    L = H.L
    D = [H.D[i] @ Minv_blocks[i] for i in range(2 ** L)]
    if L == 0:
        return HSS(0, H.rb, H.cb, D, [None], [None], [], dtype=D[0].dtype)
    Vleaf = [Minv_blocks[i].conj().T @ H.V[L][i] for i in range(2 ** L)]
    U = [None] + [list(H.U[l]) for l in range(1, L + 1)]
    V = [None] + [list(H.V[l]) for l in range(1, L)] + [Vleaf]
    B = [list(H.B[l]) for l in range(L)]
    return HSS(L, H.rb, H.cb, D, U, V, B, dtype=np.result_type(H.dtype, Minv_blocks[0].dtype))


def weighted_minnorm(H, b, w=None, Lblocks=None, slack='relax', return_factor=False):
    """argmin ||L x|| s.t. Hx = b with L = diag(w) (w nonzero) or
    L = blkdiag(Lblocks) (each invertible, conforming to the leaf columns)."""
    if w is not None:
        HL = H.scale_columns(1.0 / np.asarray(w))
        F = MinNormFactor(HL, slack=slack)
        y = F.solve(b)
        x = (y.T / np.asarray(w)).T if y.ndim == 2 else y / np.asarray(w)
    else:
        Linv = [np.linalg.inv(Lt) for Lt in Lblocks]
        HL = right_blockdiag(H, Linv)
        F = MinNormFactor(HL, slack=slack)
        y = F.solve(b)
        cb = H.cb[H.L]
        x = np.concatenate([Linv[i] @ y[cb[i]:cb[i + 1]] for i in range(2 ** H.L)], axis=0)
    return (x, F) if return_factor else x


def augment(H, lam):
    """[H, lam I] with the identity columns of leaf tau's rows placed right
    after leaf tau's own columns (a column permutation), which keeps the
    matrix HSS on the same tree:  D_tau -> [D_tau, lam I],  V_tau -> [V_tau; 0]."""
    L = H.L
    rbL, cbL = H.rb[L], H.cb[L]
    D, Vleaf, newcb = [], [], [0]
    perm_x, perm_s = [], []   # positions of x-entries and s-entries in the augmented vector
    off = 0
    for i in range(2 ** L):
        n_i, p_i = H.D[i].shape
        D.append(np.hstack([H.D[i], lam * np.eye(n_i, dtype=H.D[i].dtype)]))
        if L > 0:
            Vi = H.V[L][i]
            Vleaf.append(np.vstack([Vi, np.zeros((n_i, Vi.shape[1]), dtype=Vi.dtype)]))
        perm_x.append(np.arange(off, off + p_i))
        perm_s.append(np.arange(off + p_i, off + p_i + n_i))
        off += p_i + n_i
        newcb.append(off)
    newcb = np.array(newcb)
    cb = [None] * (L + 1)
    cb[L] = newcb
    for l in range(L - 1, -1, -1):
        cb[l] = cb[l + 1][::2]
    if L == 0:
        Ha = HSS(0, H.rb, cb, D, [None], [None], [], dtype=H.dtype)
    else:
        U = [None] + [list(H.U[l]) for l in range(1, L + 1)]
        V = [None] + [list(H.V[l]) for l in range(1, L)] + [Vleaf]
        B = [list(H.B[l]) for l in range(L)]
        Ha = HSS(L, H.rb, cb, D, U, V, B, dtype=H.dtype)
    return Ha, np.concatenate(perm_x), np.concatenate(perm_s)


def tikhonov(H, b, lam, w=None, s=None, slack='relax', return_all=False):
    """argmin ||S(Hx-b)||^2 + lam^2 ||L x||^2, L = diag(w), S = diag(s)
    (None = identity), lam > 0.  Works for wide, square and tall H."""
    Ht = H
    bt = np.asarray(b)
    if s is not None:
        Ht = Ht.scale_rows(np.asarray(s))
        bt = (bt.T * np.asarray(s)).T
    if w is not None:
        Ht = Ht.scale_columns(1.0 / np.asarray(w))
    Ha, px, ps = augment(Ht, lam)
    F = MinNormFactor(Ha, slack=slack)
    z = F.solve(bt)
    y = z[px]
    x = y if w is None else ((y.T / np.asarray(w)).T)
    if return_all:
        return x, z[ps], F, Ha
    return x
