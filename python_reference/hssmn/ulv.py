"""
Recursive ULV-type minimum-norm solver for wide HSS matrices -- an
operation-for-operation port of @hss/private/hss_ulvminnormsolve.m, plus a
factor/solve split for repeated right-hand sides.

Per leaf tau (all leaves of the current level at once):
  1. size reduction   [V^*; D] Omega = [Vt^* 0; Dt 0]           (Omega unitary)
  2. row compression  P^* U = [0; Uhat]                         (QL, P unitary)
  3. column decouple  (P^* Dt)(1:t,:) Q = [Lt 0]                (LQ, Q unitary)
  4. forced solve     Lt X1 = bhat(1:t)
  5. RHS update       bbar = (P^* b - Hhat [X1; 0])(bottom rows)   (full HSS product)
  6. discard + merge  -> reduced HSS with one level fewer; recurse
  7. reconstruct      x = Omega [Q [X1; X2]; 0]

slack='collapse' reproduces the legacy MATLAB file @hss/legacy/hss_ulvminnormsolve.m
(collapse the leaf level while some leaf violates l_tau + n_tau <= p_tau);
slack='relax' (default, as in the current MATLAB solver) keeps the tree and
uses p'_tau = min(p_tau, l_tau + n_tau) (Assumption 1 not needed).
"""
import numpy as np
from scipy.linalg import qr, solve_triangular
from .hss import HSS, _cat

EPS = np.finfo(float).eps


def ql(A):
    """A = Q L with L = [0; Lhat] (zeros on top), Q unitary (MATLAB ql())."""
    Qr, Rr = qr(A[::-1, ::-1], mode='full')
    return Qr[::-1, ::-1], Rr[::-1, ::-1]


def slack_ok(H):
    L = H.L
    if L == 0:
        return True
    for i in range(2 ** L):
        n, p = H.D[i].shape
        if H.V[L][i].shape[1] + n > p:
            return False
    return True


def _pinv_solve(D, b):
    """root solve exactly as the legacy MATLAB file: pinv for wide, backslash else."""
    m, n = D.shape
    if m < n:
        Uu, s, Vh = np.linalg.svd(D, full_matrices=False)
        tol = max(m, n) * (s[0] if s.size else 0.0) * EPS
        r = int(np.sum(s > tol))
        return Vh[:r].conj().T @ ((Uu[:, :r].conj().T @ b) / s[:r, None] if b.ndim == 2
                                  else (Uu[:, :r].conj().T @ b) / s[:r])
    return np.linalg.solve(D, b) if m == n else np.linalg.lstsq(D, b, rcond=None)[0]


class LevelFactor:
    """Factors of one recursion level."""
    __slots__ = ('L', 'rb', 'cb', 'Omega', 'P', 'Q', 'Lt', 't', 'pred',
                 'Dh1', 'Uz', 'Vh1', 'Hhat_int', 'nrows', 'ncols_red')


def factor_level(H, slack='relax'):
    """Process all leaves of H (L >= 1). Returns (LevelFactor, reduced HSS)."""
    L = H.L
    nl = 2 ** L
    F = LevelFactor()
    F.L = L
    F.Omega, F.P, F.Q, F.Lt, F.t, F.pred = [], [], [], [], [], []
    F.Dh1, F.Uz, F.Vh1 = [], [], []
    Dh_list, Uhat_list, Vh_list = [], [], []
    for i in range(nl):
        D = H.D[i]
        U = H.U[L][i]
        V = H.V[L][i]
        n, p = D.shape
        k, ell = U.shape[1], V.shape[1]
        if slack == 'strict' and ell + n > p:
            raise ValueError('insufficient slack at leaf %d' % i)
        # 1. size reduction: S = [V^*; D]; S^* = Q1 R (economic) => S Q1 = R^*
        S_star = np.hstack([V, D.conj().T])                  # p x (ell+n)
        Q1, R = qr(S_star, mode='economic')                  # Q1: p x p'
        SR = R.conj().T                                      # (ell+n) x p'
        pr = Q1.shape[1]
        Vt_star = SR[:ell, :]
        Dt = SR[ell:, :]
        # 2. row compression (QL of U)
        t = max(n - k, 0)
        if k > 0 and t > 0:
            P, LU = ql(U)
            Dt = P.conj().T @ Dt
        else:
            P, LU = None, U.copy()
        # 3. column decoupling (LQ of the top t rows)
        if t > 0:
            if t > pr:
                raise np.linalg.LinAlgError(
                    'leaf %d: t=%d > p\'=%d -- H is not of full row rank' % (i, t, pr))
            Qc, _ = qr(Dt[:t, :].conj().T, mode='full')      # p' x p'
            Dh = Dt @ Qc
            Vh_star = Vt_star @ Qc
            Lt = np.tril(Dh[:t, :t])
        else:
            Qc, Dh, Vh_star, Lt = None, Dt, Vt_star, np.zeros((0, 0))
        F.Omega.append(Q1)
        F.P.append(P)
        F.Q.append(Qc)
        F.Lt.append(Lt)
        F.t.append(t)
        F.pred.append(pr)
        # pieces needed by the RHS update  Hhat [X1; 0]
        F.Dh1.append(Dh[:, :t])
        F.Uz.append(LU)                       # n x k, top t rows zero
        F.Vh1.append(Vh_star[:, :t].conj().T)  # t x ell  (rows of Vhat for X1 columns)
        Dh_list.append(Dh)
        Uhat_list.append(LU[t:, :])
        Vh_list.append(Vh_star[:, t:].conj().T)   # (p'-t) x ell
    F.Hhat_int = (H.U, H.V, H.B)
    F.rb = H.rb
    F.cb = H.cb
    # ---- discard + merge: reduced HSS with L-1 levels
    Lp = L - 1
    Dn, Un, Vn = [], [], []
    for i in range(2 ** Lp):
        a, c = 2 * i, 2 * i + 1
        B12, B21 = H.B[Lp][i]
        ta, tc = F.t[a], F.t[c]
        Ua, Uc = Uhat_list[a], Uhat_list[c]
        Va, Vc = Vh_list[a], Vh_list[c]
        top = np.hstack([Dh_list[a][ta:, ta:], Ua @ B12 @ Vc.conj().T])
        bot = np.hstack([Uc @ B21 @ Va.conj().T, Dh_list[c][tc:, tc:]])
        Dn.append(np.vstack([top, bot]))
        if Lp >= 1:
            R, W = H.U[Lp][i], H.V[Lp][i]
            ka, la = Ua.shape[1], Va.shape[1]
            Un.append(_cat(Ua @ R[:ka], Uc @ R[ka:]))
            Vn.append(_cat(Va @ W[:la], Vc @ W[la:]))
    # boundaries of the reduced problem
    rr = [0]
    cc = [0]
    for i in range(2 ** Lp):
        rr.append(rr[-1] + Dn[i].shape[0])
        cc.append(cc[-1] + Dn[i].shape[1])
    rr, cc = np.array(rr), np.array(cc)
    rb = [None] * (Lp + 1)
    cb = [None] * (Lp + 1)
    rb[Lp], cb[Lp] = rr, cc
    for l in range(Lp - 1, -1, -1):
        rb[l], cb[l] = rb[l + 1][::2], cb[l + 1][::2]
    if Lp >= 1:
        Ured = [None] + [list(H.U[l]) for l in range(1, Lp)] + [Un]
        Vred = [None] + [list(H.V[l]) for l in range(1, Lp)] + [Vn]
    else:
        Ured, Vred = [None], [None]
    Bred = [list(H.B[l]) for l in range(Lp)]
    Hred = HSS(Lp, rb, cb, Dn, Ured, Vred, Bred, dtype=H.dtype)
    return F, Hred


def _rhs_update(F, X1):
    """y = Hhat [X1; 0] using only the X1 columns of the transformed matrix."""
    L = F.L
    U, V, B = F.Hhat_int
    nl = 2 ** L
    s = X1[0].shape[1]
    xt = [None] * (L + 1)
    xt[L] = [F.Vh1[i].conj().T @ X1[i] for i in range(nl)]
    for l in range(L - 1, 0, -1):
        xt[l] = [V[l][i].conj().T @ _cat(xt[l + 1][2 * i], xt[l + 1][2 * i + 1]) for i in range(2 ** l)]
    yt = [None] * (L + 1)
    for l in range(L):
        nxt = [None] * (2 ** (l + 1))
        for i in range(2 ** l):
            B12, B21 = B[l][i]
            ya = B12 @ xt[l + 1][2 * i + 1]
            yc = B21 @ xt[l + 1][2 * i]
            if l > 0:
                tmp = U[l][i] @ yt[l][i]
                ka = ya.shape[0]
                ya = ya + tmp[:ka]
                yc = yc + tmp[ka:]
            nxt[2 * i], nxt[2 * i + 1] = ya, yc
        yt[l + 1] = nxt
    return [F.Dh1[i] @ X1[i] + F.Uz[i] @ yt[L][i] for i in range(nl)]


class MinNormFactor:
    def __init__(self, H, slack='relax', root='factored'):
        # root='factored': apply the SVD factors to b (backward stable).
        # root='pinv': form the pseudo-inverse matrix explicitly and multiply,
        #             exactly as the legacy file's pinv(H.D)*b (not backward stable).
        self.levels = []
        self.collapses = []
        self.shape = H.shape
        self.dtype = H.dtype
        cur = H
        while True:
            nc = 0
            if slack == 'collapse':
                while cur.L > 0 and not slack_ok(cur):
                    cur = cur.merge_leaves()
                    nc += 1
            self.collapses.append(nc)
            if cur.L == 0:
                break
            F, cur = factor_level(cur, slack=slack)
            self.levels.append(F)
        self.root = cur.D[0]
        self.root_shape = cur.D[0].shape
        # precompute the root solve (current MATLAB solver; root='pinv': legacy file)
        m, n = self.root.shape
        if m < n:
            Uu, s, Vh = np.linalg.svd(self.root, full_matrices=False)
            tol = max(m, n) * (s[0] if s.size else 0.0) * EPS
            r = int(np.sum(s > tol))
            if root == 'pinv':
                P = Vh[:r].conj().T @ (Uu[:, :r].conj().T / s[:r, None])
                self._root = ('explicit', P, None, None)
            else:
                self._root = ('pinv', Uu[:, :r], s[:r], Vh[:r])
        else:
            self._root = ('solve', None, None, None)

    def _root_solve(self, b):
        kind, Uu, s, Vh = self._root
        if kind == 'pinv':
            return Vh.conj().T @ ((Uu.conj().T @ b) / s[:, None])
        if kind == 'explicit':
            return Uu @ b
        return np.linalg.solve(self.root, b)

    def solve(self, b):
        b = np.asarray(b)
        vec = b.ndim == 1
        if vec:
            b = b[:, None]
        dt = np.result_type(b.dtype, self.dtype)
        bcur = b.astype(dt, copy=False)
        stack = []
        for F in self.levels:
            nl = 2 ** F.L
            rb = F.rb[F.L]
            bh, X1, bbar = [], [], []
            for i in range(nl):
                bi = bcur[rb[i]:rb[i + 1]]
                if F.P[i] is not None:
                    bi = F.P[i].conj().T @ bi
                t = F.t[i]
                bh.append(bi)
                if t > 0:
                    X1.append(solve_triangular(F.Lt[i], bi[:t], lower=True))
                else:
                    X1.append(np.zeros((0, b.shape[1]), dtype=dt))
            y = _rhs_update(F, X1)
            for i in range(nl):
                t = F.t[i]
                bbar.append(bh[i][t:] - y[i][t:])
            stack.append(X1)
            bcur = np.concatenate(bbar, axis=0)
        x = self._root_solve(bcur)
        for F, X1 in zip(reversed(self.levels), reversed(stack)):
            nl = 2 ** F.L
            out = []
            off = 0
            for i in range(nl):
                t, pr = F.t[i], F.pred[i]
                w = pr - t
                X2 = x[off:off + w]
                off += w
                xh = _cat(X1[i], X2)
                if F.Q[i] is not None:
                    xh = F.Q[i] @ xh
                out.append(F.Omega[i] @ xh)
            x = np.concatenate(out, axis=0)
        return x[:, 0] if vec else x


def minnorm_solve(H, b, slack='relax', root='factored'):
    return MinNormFactor(H, slack=slack, root=root).solve(b)
