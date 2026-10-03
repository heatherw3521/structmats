"""
HSS constructors.

  * id_rows / id_cols         : deterministic interpolative decompositions
                                (column-pivoted QR, relative threshold on |R_jj|/|R_11|,
                                 the same rank rule as +hssutil/inter_decompv3.m but
                                 without the randomized sketch / rank cap)
  * hss_from_dense            : nested-ID construction on the MATLAB constructor's tree
                                (block rows/columns use ALL other columns/rows) -- needs A
  * hss_from_kernel           : nested-ID construction from entry evaluation + geometry,
                                near field explicit, far field through proxy points
                                (never forms A; O(N) entry evaluations)
  * lrbd_hss                  : exact generator form of the showcase's rectsampler matrix
  * random_hss                : exact random HSS with prescribed ranks
"""
import numpy as np
from scipy.linalg import qr, solve_triangular
from .hss import HSS, matlab_tree, tree_from_leaves, _cat


# ---------------------------------------------------------------- ID kernels
def _id_cols_from_qr(R, P, k):
    n = R.shape[1]
    if k == 0:
        return np.zeros((0, n), dtype=R.dtype), np.zeros(0, dtype=np.int64)
    R11 = R[:k, :k]
    R12 = R[:k, k:]
    T = solve_triangular(R11, R12) if R12.shape[1] else np.zeros((k, 0), dtype=R.dtype)
    Zp = np.hstack([np.eye(k, dtype=R.dtype), T])
    Z = np.empty_like(Zp)
    Z[:, P] = Zp
    return Z, P[:k].copy()


def id_cols(A, tol, kmax=None):
    """A ~= A[:, cols] @ Z,  Z: k x n with Z[:, cols] = I."""
    m, n = A.shape
    if m == 0 or n == 0:
        return np.zeros((0, n), dtype=A.dtype), np.zeros(0, dtype=np.int64)
    _, R, P = qr(A, mode='economic', pivoting=True)
    d = np.abs(np.diag(R))
    if d.size == 0 or d[0] == 0:
        k = 0
    else:
        k = int(np.sum(d / d[0] >= tol))   # pivoted-QR diagonal is non-increasing
    if kmax is not None:
        k = min(k, kmax)
    return _id_cols_from_qr(R, P, k)


def id_rows(A, tol, kmax=None):
    """A ~= Z @ A[rows, :],  Z: m x k with Z[rows, :] = I."""
    Zt, rows = id_cols(A.conj().T, tol, kmax)
    return Zt.conj().T, rows


# ------------------------------------------------------- generic nested ID
class KernelSpec:
    """Entry-evaluable matrix with 1-D ordered geometry (positions along a
    curve in the complex plane, monotone in the index; cyclic=True for the
    unit circle).  far_rows(I, z) must return a matrix whose columns span the
    interactions of rows I with any column located outside the proxy circle z
    (None if the far field vanishes identically); far_cols analogously."""

    def __init__(self, m, n, entries, rpos, cpos, far_rows=None, far_cols=None,
                 cyclic=False, dtype=float, name=''):
        self.m, self.n = m, n
        self.entries = entries
        self.rpos = np.asarray(rpos, dtype=complex)
        self.cpos = np.asarray(cpos, dtype=complex)
        self.far_rows = far_rows
        self.far_cols = far_cols
        self.cyclic = cyclic
        self.dtype = dtype
        self.name = name

    def dense(self):
        return self.entries(np.arange(self.m), np.arange(self.n))


def _cluster_geom(pos):
    if pos.size == 0:
        return None, 0.0
    lo = np.array([pos.real.min(), pos.imag.min()])
    hi = np.array([pos.real.max(), pos.imag.max()])
    c = complex((lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2)
    r = float(np.max(np.abs(pos - c))) if pos.size else 0.0
    return c, max(r, 1e-300)


def _near_clusters(i, nclus, geom_of, c, rad, cyclic, exclude):
    """Clusters (not in `exclude`) whose bounding disks meet |z - c| < rad,
    found by scanning outward from cluster i along the 1-D ordering.
    Returns (sorted list of near clusters, whether any far cluster exists)."""
    near = set()
    for direction in (+1, -1):
        for step in range(1, nclus):
            j = i + direction * step
            if cyclic:
                j %= nclus
            elif j < 0 or j >= nclus:
                break
            if j in exclude:
                continue
            cj, rj = geom_of(j)
            if cj is None:            # empty cluster: harmless to include, keep scanning
                near.add(j)
                continue
            if abs(cj - c) - rj < rad:
                near.add(j)
            else:
                break
    nfar = nclus - len(near) - len(exclude)
    return sorted(near), nfar > 0


def hss_build(spec, L, rb, cb, tol, mode='proxy', eta_near=3.0, eta_proxy=2.0, nproxy=96,
              kmax=None, verbose=False):
    """Nested interpolative-decomposition construction.

    mode='full'  : block rows/columns use ALL other columns/rows (the MATLAB
                   constructor's rule; O(mn) entry evaluations).
    mode='proxy' : near field explicit (full columns at the leaf level, the
                   level-(l+1) skeletons of near clusters above it), far field
                   through nproxy proxy points on a circle of radius
                   eta_proxy*r around each cluster.
    """
    m, n = spec.m, spec.n
    D = [None] * (2 ** L)
    U = [None] * (L + 1)
    V = [None] * (L + 1)
    B = [None] * L
    srow = [None] * (L + 1)   # row skeletons (global indices)
    scol = [None] * (L + 1)
    allr = np.arange(m)
    allc = np.arange(n)
    for i in range(2 ** L):
        D[i] = spec.entries(np.arange(rb[L][i], rb[L][i + 1]), np.arange(cb[L][i], cb[L][i + 1]))
    if L == 0:
        return HSS(0, rb, cb, D, [None], [None], [], dtype=spec.dtype)
    zq = np.exp(2j * np.pi * (np.arange(nproxy) + 0.5) / nproxy)
    for l in range(L, 0, -1):
        nc = 2 ** l
        Ul, Vl, sr, sc = [None] * nc, [None] * nc, [None] * nc, [None] * nc
        rgeo = [_cluster_geom(spec.rpos[rb[l][i]:rb[l][i + 1]]) for i in range(nc)]
        cgeo = [_cluster_geom(spec.cpos[cb[l][i]:cb[l][i + 1]]) for i in range(nc)]
        if l < L:
            # geometry of level l+1 clusters (for skeleton-based near fields)
            rgeo1 = [_cluster_geom(spec.rpos[rb[l + 1][j]:rb[l + 1][j + 1]]) for j in range(2 * nc)]
            cgeo1 = [_cluster_geom(spec.cpos[cb[l + 1][j]:cb[l + 1][j + 1]]) for j in range(2 * nc)]
        for i in range(nc):
            # ---------------- candidate rows of this node
            if l == L:
                cand_r = np.arange(rb[l][i], rb[l][i + 1])
                cand_c = np.arange(cb[l][i], cb[l][i + 1])
            else:
                cand_r = np.concatenate([srow[l + 1][2 * i], srow[l + 1][2 * i + 1]])
                cand_c = np.concatenate([scol[l + 1][2 * i], scol[l + 1][2 * i + 1]])
            # ---------------- row ID: interactions with columns outside J_i
            if len(cand_r) == 0:
                Ul[i] = np.zeros((0, 0), dtype=spec.dtype)
                sr[i] = cand_r
            elif mode == 'full':
                outside = np.concatenate([allc[:cb[l][i]], allc[cb[l][i + 1]:]])
                blk = spec.entries(cand_r, outside)
            else:
                c0, r0 = rgeo[i]
                if c0 is None:
                    c0, r0 = cgeo[i]
                if l == L:
                    near, has_far = _near_clusters(i, nc, lambda j: cgeo[j], c0, eta_near * r0,
                                                   spec.cyclic, {i})
                    nearcols = np.concatenate([np.arange(cb[l][j], cb[l][j + 1]) for j in near]) \
                        if near else np.zeros(0, dtype=np.int64)
                else:
                    near, has_far = _near_clusters(2 * i, 2 * nc, lambda j: cgeo1[j], c0, eta_near * r0,
                                                   spec.cyclic, {2 * i, 2 * i + 1})
                    nearcols = np.concatenate([scol[l + 1][j] for j in near]) \
                        if near else np.zeros(0, dtype=np.int64)
                parts = [spec.entries(cand_r, nearcols)]
                if has_far and spec.far_rows is not None:
                    pr_ = spec.far_rows(cand_r, c0 + eta_proxy * r0 * zq)
                    if not np.iscomplexobj(np.zeros(0, dtype=spec.dtype)):
                        pr_ = np.hstack([pr_.real, pr_.imag])   # real kernel: real span
                    parts.append(pr_)
                blk = np.hstack(parts)
            if len(cand_r) > 0:
                Z, sel = id_rows(blk, tol, kmax)
                Ul[i] = Z
                sr[i] = cand_r[sel]
            # ---------------- column ID: interactions with rows outside I_i
            if mode == 'full':
                outside = np.concatenate([allr[:rb[l][i]], allr[rb[l][i + 1]:]])
                blk = spec.entries(outside, cand_c)
            else:
                c0, r0 = cgeo[i]
                if l == L:
                    near, has_far = _near_clusters(i, nc, lambda j: rgeo[j], c0, eta_near * r0,
                                                   spec.cyclic, {i})
                    nearrows = np.concatenate([np.arange(rb[l][j], rb[l][j + 1]) for j in near]) \
                        if near else np.zeros(0, dtype=np.int64)
                else:
                    near, has_far = _near_clusters(2 * i, 2 * nc, lambda j: rgeo1[j], c0, eta_near * r0,
                                                   spec.cyclic, {2 * i, 2 * i + 1})
                    nearrows = np.concatenate([srow[l + 1][j] for j in near]) \
                        if near else np.zeros(0, dtype=np.int64)
                parts = [spec.entries(nearrows, cand_c)]
                if has_far and spec.far_cols is not None:
                    pc_ = spec.far_cols(c0 + eta_proxy * r0 * zq, cand_c)
                    if not np.iscomplexobj(np.zeros(0, dtype=spec.dtype)):
                        pc_ = np.vstack([pc_.real, pc_.imag])
                    parts.append(pc_)
                blk = np.vstack(parts)
            Y, sel = id_cols(blk, tol, kmax)       # blk ~= blk[:, sel] Y
            Vl[i] = Y.conj().T                      # V^* = Y
            sc[i] = cand_c[sel]
        U[l], V[l], srow[l], scol[l] = Ul, Vl, sr, sc
        if verbose:
            print('level %d: max k=%d max l=%d' % (l, max(u.shape[1] for u in Ul), max(v.shape[1] for v in Vl)))
    for l in range(L):
        Bl = []
        for i in range(2 ** l):
            a, c = 2 * i, 2 * i + 1
            B12 = spec.entries(srow[l + 1][a], scol[l + 1][c])
            B21 = spec.entries(srow[l + 1][c], scol[l + 1][a])
            Bl.append((B12, B21))
        B[l] = Bl
    U[0] = None
    V[0] = None
    H = HSS(L, rb, cb, D, U, V, B, dtype=spec.dtype)
    H.skel = (srow, scol)
    return H


def hss_from_dense(A, blocksize, tol=1e-12):
    """Mirror of hss(A, blocksize=..., tol=...) (deterministic ID)."""
    m, n = A.shape
    L, rb, cb = matlab_tree(m, n, blocksize)
    spec = KernelSpec(m, n, lambda I, J: A[np.ix_(I, J)], np.arange(m), np.arange(n),
                      dtype=A.dtype)
    return hss_build(spec, L, rb, cb, tol, mode='full')


# ----------------------------------------------------- generator families
def lrbd_hss(m, n, k, blocksize, rng, tree=None):
    """F1: A = X Y + blkdiag(E_tau), X (m x k), Y (k x n), E_tau iid U[0,1];
    exactly the matrix the showcase's rectsampler() produces, but assembled
    directly in generator form (never dense): U_tau = X(I_tau,:),
    V_tau = Y(:,J_tau)^T, R = W = [I;I], B = I."""
    L, rb, cb = matlab_tree(m, n, blocksize) if tree is None else tree
    X = rng.random((m, k))
    Y = rng.random((k, n))
    I = np.eye(k)
    D = []
    for i in range(2 ** L):
        r0, r1, c0, c1 = rb[L][i], rb[L][i + 1], cb[L][i], cb[L][i + 1]
        D.append(X[r0:r1] @ Y[:, c0:c1] + rng.random((r1 - r0, c1 - c0)))
    if L == 0:
        return HSS(0, rb, cb, D, [None], [None], [])
    U = [None] * (L + 1)
    V = [None] * (L + 1)
    U[L] = [X[rb[L][i]:rb[L][i + 1]].copy() for i in range(2 ** L)]
    V[L] = [Y[:, cb[L][i]:cb[L][i + 1]].T.copy() for i in range(2 ** L)]
    for l in range(1, L):
        U[l] = [np.vstack([I, I]) for _ in range(2 ** l)]
        V[l] = [np.vstack([I, I]) for _ in range(2 ** l)]
    B = [[(I.copy(), I.copy()) for _ in range(2 ** l)] for l in range(L)]
    return HSS(L, rb, cb, D, U, V, B)


def _orth(rng, a, b, cplx=False):
    G = rng.standard_normal((a, b))
    if cplx:
        G = G + 1j * rng.standard_normal((a, b))
    Q, _ = np.linalg.qr(G)
    return Q


def random_hss(m, n, k, blocksize, rng, coupling=1.0, cplx=False, tree=None, ell=None):
    """F2: exact random HSS with ranks k (rows) / ell (cols) on every node.
    D_tau iid N(0,1)/sqrt(p_tau) (n_tau x p_tau); U, V, R, W with orthonormal
    columns; B iid N(0,1)*coupling/sqrt(k)."""
    ell = k if ell is None else ell
    L, rb, cb = matlab_tree(m, n, blocksize) if tree is None else tree
    def g(a, b):
        G = rng.standard_normal((a, b))
        return G + 1j * rng.standard_normal((a, b)) if cplx else G
    D = []
    for i in range(2 ** L):
        nr, nc = rb[L][i + 1] - rb[L][i], cb[L][i + 1] - cb[L][i]
        D.append(g(nr, nc) / np.sqrt(nc))
    if L == 0:
        return HSS(0, rb, cb, D, [None], [None], [])
    U = [None] * (L + 1)
    V = [None] * (L + 1)
    U[L] = [_orth(rng, rb[L][i + 1] - rb[L][i], k, cplx) for i in range(2 ** L)]
    V[L] = [_orth(rng, cb[L][i + 1] - cb[L][i], ell, cplx) for i in range(2 ** L)]
    for l in range(1, L):
        U[l] = [_orth(rng, 2 * k, k, cplx) for _ in range(2 ** l)]
        V[l] = [_orth(rng, 2 * ell, ell, cplx) for _ in range(2 ** l)]
    B = [[(g(k, ell) * coupling / np.sqrt(k), g(k, ell) * coupling / np.sqrt(k))
          for _ in range(2 ** l)] for l in range(L)]
    return HSS(L, rb, cb, D, U, V, B, dtype=complex if cplx else float)
