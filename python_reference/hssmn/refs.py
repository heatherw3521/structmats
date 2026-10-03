"""Reference solvers and error metrics."""
import numpy as np
from scipy.linalg import qr, solve_triangular, cho_factor, cho_solve
from scipy.sparse.linalg import LinearOperator, cg


def minnorm_qr(A, b):
    """Backward-stable dense min-norm solution via QR of A^*: A^* = Q R,
    x = Q R^{-*} b (A full row rank)."""
    Q, R = qr(A.conj().T, mode='economic')
    return Q @ solve_triangular(R, b, trans='C', lower=False)


def minnorm_pinv(A, b):
    return np.linalg.pinv(A) @ b


def minnorm_normal_eq(A, b):
    """x = A^* (A A^*)^{-1} b via Cholesky (squares the condition number)."""
    G = A @ A.conj().T
    c = cho_factor(G)
    return A.conj().T @ cho_solve(c, b)


def minnorm_cgne(H, b, rtol=1e-14, maxiter=None):
    """CG on H H^* y = b (Craig/CGNE) using HSS matvecs; x = H^* y."""
    m, n = H.shape
    Ha = H.adjoint()
    dt = np.result_type(H.dtype, b.dtype)
    op = LinearOperator((m, m), matvec=lambda y: H.matvec(Ha.matvec(y)), dtype=dt)
    it = [0]
    def cb(_):
        it[0] += 1
    y, info = cg(op, b, rtol=rtol, maxiter=maxiter or 20 * m, callback=cb)
    return Ha.matvec(y), it[0], info


def minnorm_mp(A, b, dps=50):
    """High-precision reference x = A^T (A A^T)^{-1} b (mpmath)."""
    import mpmath as mp
    mp.mp.dps = dps
    Am = mp.matrix(A.tolist())
    bm = mp.matrix(b.reshape(-1, 1).tolist())
    At = Am.T if not np.iscomplexobj(A) else Am.H
    G = Am * At
    y = mp.lu_solve(G, bm)
    x = At * y
    return np.array([complex(v) if np.iscomplexobj(A) else float(v) for v in x])


def tikhonov_svd(A, b, lam):
    U, s, Vh = np.linalg.svd(A, full_matrices=False)
    return Vh.conj().T @ ((s / (s ** 2 + lam ** 2)) * (U.conj().T @ b))


def rel(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def backward_error(H_matvec, normH, x, b):
    """normwise (Rigal-Gaches) backward error ||b - Hx||/(||H|| ||x|| + ||b||)."""
    r = b - H_matvec(x)
    return float(np.linalg.norm(r) / (normH * np.linalg.norm(x) + np.linalg.norm(b)))


def nullspace_fraction(A, x):
    """||P_null(A) x|| / ||x|| computed densely (QR of A^*)."""
    Q, _ = qr(A.conj().T, mode='economic')
    xr = Q @ (Q.conj().T @ x)
    return float(np.linalg.norm(x - xr) / np.linalg.norm(x))


def hss_norm2_est(H, iters=30, rng=None):
    """power iteration estimate of ||H||_2 with HSS matvecs"""
    rng = np.random.default_rng(0) if rng is None else rng
    m, n = H.shape
    Ha = H.adjoint()
    v = rng.standard_normal(n)
    if np.iscomplexobj(np.zeros(0, dtype=H.dtype)):
        v = v + 1j * rng.standard_normal(n)
    v /= np.linalg.norm(v)
    s = 0.0
    for _ in range(iters):
        w = Ha.matvec(H.matvec(v))
        s = np.linalg.norm(w)
        v = w / s
    return float(np.sqrt(s))
