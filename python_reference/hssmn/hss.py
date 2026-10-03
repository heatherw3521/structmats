"""
Generator-form HSS matrices on a perfect binary tree (reference implementation).

Conventions (match the memo; code correspondence to the MATLAB @hss class in
brackets):

  * levels l = 0 (root) ... L (leaves); node (l, i), i = 0 .. 2^l - 1;
    children of (l, i) are (l+1, 2i) and (l+1, 2i+1).
  * rb[l], cb[l]: row / column boundaries of the level-l clusters
    (length 2^l + 1, global offsets).
  * leaves:      D[i]      n_i x p_i          [leaf .D]
                 U[L][i]   n_i x k_i          [Z of the off-diagonal block in the leaf's row]
                 V[L][i]   p_i x l_i          [Y' of the off-diagonal block in the leaf's column]
  * internal non-root nodes:
                 U[l][i] = R_i  (k_{2i}+k_{2i+1}) x k_i    (translation, nested row basis)
                 V[l][i] = W_i  (l_{2i}+l_{2i+1}) x l_i    (translation, nested column basis)
  * couplings stored at the parent (l < L):
                 B[l][i] = (B12, B21),  B12: k_{l+1,2i} x l_{l+1,2i+1}
                                        B21: k_{l+1,2i+1} x l_{l+1,2i}       [lrcomponent]
    so that  H(I_a, J_c) = Ubig_a B12 Vbig_c^*  for the children a=2i, c=2i+1.
"""
import numpy as np


def _cat(a, b):
    return np.concatenate((a, b), axis=0)


class HSS:
    def __init__(self, L, rb, cb, D, U, V, B, dtype=float):
        self.L = L
        self.rb = [np.asarray(r, dtype=np.int64) for r in rb]
        self.cb = [np.asarray(c, dtype=np.int64) for c in cb]
        self.D = D
        self.U = U
        self.V = V
        self.B = B
        self.dtype = np.result_type(*([d.dtype for d in D] + [dtype]))

    # ------------------------------------------------------------------ sizes
    @property
    def shape(self):
        return (int(self.rb[0][-1] - self.rb[0][0]), int(self.cb[0][-1] - self.cb[0][0]))

    def leaf_sizes(self):
        L = self.L
        return [(d.shape[0], d.shape[1]) for d in self.D]

    def rank_k(self, l, i):
        return self.U[l][i].shape[1]

    def rank_l(self, l, i):
        return self.V[l][i].shape[1]

    def max_rank(self):
        r = 0
        for l in range(1, self.L + 1):
            for i in range(2 ** l):
                r = max(r, self.U[l][i].shape[1], self.V[l][i].shape[1])
        return r

    def storage(self):
        """number of stored scalars"""
        s = sum(d.size for d in self.D)
        for l in range(1, self.L + 1):
            for i in range(2 ** l):
                s += self.U[l][i].size + self.V[l][i].size
        for l in range(self.L):
            for i in range(2 ** l):
                s += self.B[l][i][0].size + self.B[l][i][1].size
        return s

    # ---------------------------------------------------------------- matvec
    def matvec(self, x):
        """y = H x  (x: n or n x s)"""
        x = np.asarray(x)
        vec = x.ndim == 1
        if vec:
            x = x[:, None]
        L = self.L
        if L == 0:
            y = self.D[0] @ x
            return y[:, 0] if vec else y
        cb, rb = self.cb[L], self.rb[L]
        nl = 2 ** L
        dt = np.result_type(self.dtype, x.dtype)
        # upward pass: xt[l][i] = Vbig^* x_i
        xt = [None] * (L + 1)
        xt[L] = [self.V[L][i].conj().T @ x[cb[i]:cb[i + 1]] for i in range(nl)]
        for l in range(L - 1, 0, -1):
            xt[l] = [self.V[l][i].conj().T @ _cat(xt[l + 1][2 * i], xt[l + 1][2 * i + 1])
                     for i in range(2 ** l)]
        # downward pass: yt[l][i] coefficient in Ubig_i's basis
        yt = [None] * (L + 1)
        for l in range(L):
            nxt = [None] * (2 ** (l + 1))
            for i in range(2 ** l):
                B12, B21 = self.B[l][i]
                ya = B12 @ xt[l + 1][2 * i + 1]
                yc = B21 @ xt[l + 1][2 * i]
                if l > 0:
                    tmp = self.U[l][i] @ yt[l][i]
                    ka = ya.shape[0]
                    ya = ya + tmp[:ka]
                    yc = yc + tmp[ka:]
                nxt[2 * i] = ya
                nxt[2 * i + 1] = yc
            yt[l + 1] = nxt
        y = np.empty((self.shape[0], x.shape[1]), dtype=dt)
        for i in range(nl):
            y[rb[i]:rb[i + 1]] = self.D[i] @ x[cb[i]:cb[i + 1]] + self.U[L][i] @ yt[L][i]
        return y[:, 0] if vec else y

    def __matmul__(self, x):
        return self.matvec(x)

    def adjoint(self):
        """H^* as an HSS matrix (row and column roles swapped)."""
        L = self.L
        D = [d.conj().T for d in self.D]
        U = [None] + [list(self.V[l]) for l in range(1, L + 1)]
        V = [None] + [list(self.U[l]) for l in range(1, L + 1)]
        B = [[(b[1].conj().T, b[0].conj().T) for b in self.B[l]] for l in range(L)]
        return HSS(L, self.cb, self.rb, D, U, V, B, dtype=self.dtype)

    def rmatvec(self, y):
        """H^* y"""
        return self.adjoint().matvec(y)

    # ----------------------------------------------------------------- dense
    def _ubig(self, l, i):
        if l == self.L:
            return self.U[l][i]
        Ua = self._ubig(l + 1, 2 * i)
        Uc = self._ubig(l + 1, 2 * i + 1)
        R = self.U[l][i]
        ka = Ua.shape[1]
        return _cat(Ua @ R[:ka], Uc @ R[ka:])

    def _vbig(self, l, i):
        if l == self.L:
            return self.V[l][i]
        Va = self._vbig(l + 1, 2 * i)
        Vc = self._vbig(l + 1, 2 * i + 1)
        W = self.V[l][i]
        la = Va.shape[1]
        return _cat(Va @ W[:la], Vc @ W[la:])

    def dense(self):
        m, n = self.shape
        A = np.zeros((m, n), dtype=self.dtype)
        L = self.L
        for i in range(2 ** L):
            A[self.rb[L][i]:self.rb[L][i + 1], self.cb[L][i]:self.cb[L][i + 1]] = self.D[i]
        for l in range(L):
            for i in range(2 ** l):
                a, c = 2 * i, 2 * i + 1
                B12, B21 = self.B[l][i]
                ra = slice(self.rb[l + 1][a], self.rb[l + 1][a + 1])
                rc = slice(self.rb[l + 1][c], self.rb[l + 1][c + 1])
                ca = slice(self.cb[l + 1][a], self.cb[l + 1][a + 1])
                cc = slice(self.cb[l + 1][c], self.cb[l + 1][c + 1])
                A[ra, cc] = self._ubig(l + 1, a) @ B12 @ self._vbig(l + 1, c).conj().T
                A[rc, ca] = self._ubig(l + 1, c) @ B21 @ self._vbig(l + 1, a).conj().T
        return A

    # ----------------------------------------------------- structural edits
    def merge_leaves(self):
        """Collapse the leaf level into the parents (MATLAB `levelup`): the
        represented matrix is unchanged, the tree loses one level."""
        L = self.L
        assert L >= 1
        D, Un, Vn = [], [], []
        for i in range(2 ** (L - 1)):
            a, c = 2 * i, 2 * i + 1
            B12, B21 = self.B[L - 1][i]
            Ua, Uc = self.U[L][a], self.U[L][c]
            Va, Vc = self.V[L][a], self.V[L][c]
            top = np.hstack([self.D[a], Ua @ B12 @ Vc.conj().T])
            bot = np.hstack([Uc @ B21 @ Va.conj().T, self.D[c]])
            D.append(np.vstack([top, bot]))
            if L - 1 >= 1:
                R, W = self.U[L - 1][i], self.V[L - 1][i]
                ka, la = Ua.shape[1], Va.shape[1]
                Un.append(_cat(Ua @ R[:ka], Uc @ R[ka:]))
                Vn.append(_cat(Va @ W[:la], Vc @ W[la:]))
        U = [list(self.U[l]) if l > 0 else None for l in range(L - 1)] + ([Un] if L - 1 >= 1 else [])
        V = [list(self.V[l]) if l > 0 else None for l in range(L - 1)] + ([Vn] if L - 1 >= 1 else [])
        if L - 1 == 0:
            U, V = [None], [None]
        B = [list(self.B[l]) for l in range(L - 1)]
        return HSS(L - 1, self.rb[:L], self.cb[:L], D, U, V, B, dtype=self.dtype)

    def scale_columns(self, s):
        """H diag(s) for a column scaling s (vector) -- stays HSS on the same tree."""
        L = self.L
        cb = self.cb[L]
        D = [self.D[i] * s[cb[i]:cb[i + 1]][None, :] for i in range(2 ** L)]
        if L == 0:
            return HSS(0, self.rb, self.cb, D, [None], [None], [], dtype=self.dtype)
        V = [None] + [list(self.V[l]) for l in range(1, L)] + \
            [[np.conj(s[cb[i]:cb[i + 1]])[:, None] * self.V[L][i] for i in range(2 ** L)]]
        U = [None] + [list(self.U[l]) for l in range(1, L + 1)]
        B = [list(self.B[l]) for l in range(L)]
        return HSS(L, self.rb, self.cb, D, U, V, B, dtype=np.result_type(self.dtype, s.dtype))

    def scale_rows(self, s):
        """diag(s) H"""
        return self.adjoint().scale_columns(np.conj(s)).adjoint()

    def copy(self):
        L = self.L
        return HSS(L, self.rb, self.cb, [d.copy() for d in self.D],
                   [None] + [[u.copy() for u in self.U[l]] for l in range(1, L + 1)],
                   [None] + [[v.copy() for v in self.V[l]] for l in range(1, L + 1)],
                   [[(b[0].copy(), b[1].copy()) for b in self.B[l]] for l in range(L)],
                   dtype=self.dtype)

    def summary(self):
        m, n = self.shape
        ns = [d.shape for d in self.D]
        return dict(m=m, n=n, L=self.L, leaf_rows=(min(s[0] for s in ns), max(s[0] for s in ns)),
                    leaf_cols=(min(s[1] for s in ns), max(s[1] for s in ns)),
                    max_rank=self.max_rank() if self.L > 0 else 0, storage=self.storage())


# ---------------------------------------------------------------- partitions
def matlab_tree(m, n, blocksize):
    """Tree exactly as @hss/private/hss_constructor.m builds it:
    depth = largest level at which the worst-case (floor) branch still has
    both dimensions >= blocksize; first child gets ceil(size/2)."""
    L = 0
    mm, nn = m, n
    while True:
        m2, n2 = mm // 2, nn // 2
        if m2 < blocksize or n2 < blocksize:
            break
        mm, nn = m2, n2
        L += 1
    rb = [np.array([0, m])]
    cb = [np.array([0, n])]
    for l in range(L):
        r, c = [0], [0]
        for i in range(2 ** l):
            r0, r1 = rb[l][i], rb[l][i + 1]
            c0, c1 = cb[l][i], cb[l][i + 1]
            r += [r0 + int(np.ceil((r1 - r0) / 2)), r1]
            c += [c0 + int(np.ceil((c1 - c0) / 2)), c1]
        rb.append(np.array(r))
        cb.append(np.array(c))
    return L, rb, cb


def tree_from_leaves(row_leaf_bounds, col_leaf_bounds):
    """Given leaf-level boundaries (length 2^L + 1), build all levels."""
    rbL = np.asarray(row_leaf_bounds)
    cbL = np.asarray(col_leaf_bounds)
    nl = len(rbL) - 1
    L = int(round(np.log2(nl)))
    assert 2 ** L == nl
    rb = [None] * (L + 1)
    cb = [None] * (L + 1)
    rb[L], cb[L] = rbL, cbL
    for l in range(L - 1, -1, -1):
        rb[l] = rb[l + 1][::2]
        cb[l] = cb[l + 1][::2]
    return L, rb, cb
