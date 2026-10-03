"""Flop count of one factorization (memo 1, Sec. 5) for an F2-type matrix:
ranks k (rows) = l (cols) at every node, perfect tree of depth L, relaxed width.
Counts the dense kernels the factorization performs per node: economic QR of the
size-reduction stack and its Q, QL of U with the full P, P^* D, LQ of the top rows
with the full Q, and the products D Q, V^* Q.  LAPACK flop counts."""


def _qr(m, n):
    n = min(m, n)
    return 2.0 * n * n * (m - n / 3.0)


def _qfull(m, k):
    return 4.0 * m * m * k - 4.0 * m * k * k + 4.0 / 3.0 * k ** 3


def _qecon(m, k):
    return 2.0 * m * k * k - 2.0 / 3.0 * k ** 3


def node_flops(nr, nc, k, l):
    c = min(nc, l + nr)                  # size-reduced width p'
    t = max(nr - k, 0)                   # decoupled rows
    f = _qr(nc, c) + _qecon(nc, c)       # size reduction
    f += _qr(nr, k) + _qfull(nr, k) + 2.0 * nr * nr * c      # QL of U, P^* D
    f += _qr(c, t) + _qfull(c, t) + 2.0 * nr * c * c + 2.0 * l * c * c   # decoupling
    return f, k, c - t


def factor_flops(m, n, L, k, l=None):
    l = k if l is None else l
    nr, nc, nodes = m // 2 ** L, n // 2 ** L, 2 ** L
    F = 0.0
    for _ in range(L):
        f, kr, fc = node_flops(nr, nc, k, l)
        F += nodes * f
        nodes //= 2
        nr, nc = 2 * kr, 2 * fc
    F += 6.0 * nr * nr * nc              # root SVD (order of magnitude)
    return F
