"""
Test-matrix families used in the memos.  Every family is defined by explicit
formulas so that the exact matrix can be regenerated in MATLAB (see the
`hss_examples/minnorm_paper` suite) -- see the memo's test-matrix table.
"""
import numpy as np
from .build import KernelSpec, hss_build, lrbd_hss, random_hss
from .hss import matlab_tree


# ----------------------------------------------------------------- F3 Cauchy
def cauchy_interlaced_points(M, ratio=3):
    """x_i = i (i=1..M); `ratio` points y strictly inside each gap
    (y = x_i + q/(ratio+1), q=1..ratio)  ->  N = (M-1)*ratio.
    Identical to cauchyInterlaced() in hss_minnorm_showcase.m."""
    x = np.arange(1, M + 1, dtype=float)
    q = np.arange(1, ratio + 1) / (ratio + 1)
    y = (x[:-1, None] + q[None, :]).ravel()
    return x, y


def cauchy_spec(M, ratio=3):
    x, y = cauchy_interlaced_points(M, ratio)
    ent = lambda I, J: 1.0 / (x[I][:, None] - y[J][None, :])
    far_r = lambda I, z: 1.0 / (x[I][:, None] - z[None, :])
    far_c = lambda z, J: 1.0 / (z[:, None] - y[J][None, :])
    return KernelSpec(len(x), len(y), ent, x, y, far_r, far_c, cyclic=False, dtype=float,
                      name='cauchy_interlaced')


# --------------------------------------------- F4 decimated convolution
def wavelet(kind, sigma, w):
    d = np.arange(-w, w + 1, dtype=float)
    if kind == 'gauss':
        g = np.exp(-d ** 2 / (2 * sigma ** 2))
    elif kind == 'ricker':
        g = (1 - (d / sigma) ** 2) * np.exp(-d ** 2 / (2 * sigma ** 2))
    else:
        raise ValueError(kind)
    return g


def conv_spec(n, stride=2, kind='gauss', sigma=2.0, halfwidth=None, mode='valid'):
    """Decimated ('strided') convolution  y_i = sum_j g(s*i + o - j) x_j.
    g is truncated to |d| <= w (exactly banded).  mode='valid': only outputs
    whose stencil lies inside [0, n-1] (o = w);  m = floor((n-1-2w)/s)+1."""
    w = int(np.ceil(4 * sigma)) if halfwidth is None else halfwidth
    g = wavelet(kind, sigma, w)
    if mode == 'valid':
        o = w
        m = (n - 1 - 2 * w) // stride + 1
    else:  # 'same' (zero boundary)
        o = 0
        m = (n - 1) // stride + 1
    rpos = (stride * np.arange(m) + o).astype(float)
    cpos = np.arange(n, dtype=float)

    def ent(I, J):
        d = (rpos[I][:, None] - cpos[J][None, :]).astype(np.int64)
        out = np.zeros(d.shape)
        mask = np.abs(d) <= w
        out[mask] = g[d[mask] + w]
        return out
    spec = KernelSpec(m, n, ent, rpos, cpos, None, None, cyclic=False, dtype=float,
                      name='conv_%s' % kind)
    spec.w, spec.g, spec.stride, spec.offset = w, g, stride, o
    return spec


def conv_apply(spec, x):
    """exact y = A x via np.convolve (for verification at any size)."""
    full = np.convolve(x, spec.g[::-1], mode='full')   # full[k] = sum_j g(k - w - j)... see below
    # y_i = sum_j g(r_i - j) x_j with r_i = s*i + o ; np.convolve(x, g_rev)[k] = sum_j x_j g_rev[k-j]
    # g_rev[t] = g(w - t) -> term g(w - k + j); need w - k + j = -(r_i - j) ... use direct form:
    w = spec.w
    c = np.convolve(x, spec.g, mode='full')            # c[k] = sum_j x_j g[k-j] = sum_j x_j g(k-j-w)
    r = (spec.stride * np.arange(spec.m) + spec.offset)
    return c[r + w]


# --------------------------------------------------- F5 NUDFT Cauchy-like
def nudft_points(m, rng, jitter=0.5):
    """jittered points in [0,1):  x_j = (j + 1/2 + jitter*u_j)/m, u ~ U[-1/2,1/2]."""
    j = np.arange(m)
    return np.sort((j + 0.5 + jitter * (rng.random(m) - 0.5)) / m)


def nudft_spec(xs, N):
    """C = A F^*,  A_jk = exp(2 pi i k x_j) (k = -N/2..N/2-1), F unitary DFT on
    the grid l/N.  C_jl = D_N(x_j - l/N)/sqrt(N),
    D_N(u) = e^{-i pi u} sin(pi N u)/sin(pi u).  Cauchy-like:
    C_jl = a_j b_l /(z_j - w_l),  z = e^{2 pi i x}, w = e^{2 pi i l/N}."""
    m = len(xs)
    grid = np.arange(N) / N
    z = np.exp(2j * np.pi * xs)
    wl = np.exp(2j * np.pi * grid)
    # a_j b_l/(z_j - w_l) with sin(pi(x-y)) = (z - w) e^{-i pi (x+y)}/(2i)
    # => D_N(u)/sqrt(N) = e^{-i pi u} sin(pi N u) 2i e^{i pi (x+y)}/((z-w) sqrt(N))
    # sin(pi N (x - l/N)) = (-1)^l sin(pi N x)
    sgn = (-1.0) ** np.arange(N)
    a = 2j * np.sin(np.pi * N * xs) / np.sqrt(N)          # times e^{-i pi u} e^{i pi(x+y)} = e^{2 i pi y}
    b = sgn * np.exp(2j * np.pi * grid)                    # e^{2 i pi y} (-1)^l
    def ent(I, J):
        u = xs[I][:, None] - grid[J][None, :]
        num = np.sin(np.pi * N * u)
        den = np.sin(np.pi * u)
        out = np.exp(-1j * np.pi * u) * num / (den * np.sqrt(N))
        small = np.abs(den) < 1e-14
        if np.any(small):
            out[small] = np.sqrt(N)
        return out
    far_r = lambda I, zq: a[I][:, None] / (z[I][:, None] - zq[None, :])
    far_c = lambda zq, J: b[J][None, :] / (zq[:, None] - wl[J][None, :])
    spec = KernelSpec(m, N, ent, z, wl, far_r, far_c, cyclic=True, dtype=complex,
                      name='nudft_cauchy')
    spec.xs, spec.N, spec.a, spec.b, spec.z, spec.w = xs, N, a, b, z, wl
    return spec


def nudft_vandermonde(xs, N):
    k = np.arange(-N // 2, N // 2)
    return np.exp(2j * np.pi * np.outer(xs, k))


def unitary_dft(N):
    """F with (F c)_l = sum_k c_k e^{2 pi i k l/N}/sqrt(N), k = -N/2..N/2-1,
    so that A F^* = C  (A = NUDFT Vandermonde)."""
    k = np.arange(-N // 2, N // 2)
    l = np.arange(N)
    return np.exp(2j * np.pi * np.outer(l, k) / N) / np.sqrt(N)


# ------------------------------------- F6 Toeplitz -> Cauchy-like (rank 2)
class ToeplitzCauchy:
    """T (m x n) Toeplitz, T_ij = t_{i-j};  tc = T[:,0], tr = T[0,:].
    With Z_theta = cyclic down-shift with theta in the corner,
        Z_1^{(m)} T - T Z_theta^{(n)} = e_1 r^T + c e_n^T          (rank <= 2)
    and the unitary transforms F_m (F_jk = e^{-2 pi i jk/m}/sqrt m),
    D = diag(theta^{k/n}),
        C = F_m T D^* F_n^*,   Lam C - C Mu = G H^*,
        lam_j = e^{-2 pi i j/m},  mu_l = theta^{1/n} e^{-2 pi i l/n},
        G = F_m [e_1, c],  H = F_n D [conj(r), e_n].
    T x = b  <=>  C y = F_m b,  y = F_n D x  (norm preserving)."""

    def __init__(self, tc, tr, phi=np.pi):
        tc = np.asarray(tc, dtype=complex)
        tr = np.asarray(tr, dtype=complex)
        assert tc[0] == tr[0]
        m, n = len(tc), len(tr)
        self.m, self.n, self.phi = m, n, phi
        self.tc, self.tr = tc, tr
        theta = np.exp(1j * phi)
        def T_entry(i, j):
            d = i - j
            return tc[d] if d >= 0 else tr[-d]
        # row 1 of the displacement: T_{m,j} - T_{1,j+1} (j<n), T_{m,n} - theta T_{1,1}
        lastrow = np.array([T_entry(m - 1, jj) for jj in range(n)])
        r = lastrow - np.concatenate([tr[1:], [theta * tr[0]]])
        # column n, rows 2..m: T_{i-1,n} - theta T_{i,1}
        lastcol = np.array([T_entry(ii, n - 1) for ii in range(m)])
        c = np.zeros(m, dtype=complex)
        c[1:] = lastcol[:-1] - theta * tc[1:]
        self.r, self.c = r, c
        self.d = np.exp(1j * phi * np.arange(n) / n)          # D = diag(theta^{k/n})
        self.lam = np.exp(-2j * np.pi * np.arange(m) / m)
        self.mu = np.exp(1j * phi / n) * np.exp(-2j * np.pi * np.arange(n) / n)
        e1 = np.zeros(m, dtype=complex); e1[0] = 1
        en = np.zeros(n, dtype=complex); en[-1] = 1
        self.G = np.stack([self.Fm(e1), self.Fm(c)], axis=1)
        self.H = np.stack([self.FnD(np.conj(r)), self.FnD(en)], axis=1)

    # unitary transforms via FFT (numpy fft uses e^{-2 pi i jk/N})
    def Fm(self, v):
        return np.fft.fft(v, axis=0) / np.sqrt(self.m)

    def Fm_inv(self, v):
        return np.fft.ifft(v, axis=0) * np.sqrt(self.m)

    def FnD(self, x):
        dd = self.d if x.ndim == 1 else self.d[:, None]
        return np.fft.fft(dd * x, axis=0) / np.sqrt(self.n)

    def FnD_inv(self, y):
        dd = self.d if y.ndim == 1 else self.d[:, None]
        return np.conj(dd) * np.fft.ifft(y, axis=0) * np.sqrt(self.n)

    def T_dense(self):
        from scipy.linalg import toeplitz
        return toeplitz(self.tc, self.tr)

    def C_dense_fft(self):
        T = self.T_dense()
        A = self.Fm(T)                         # F_m T
        # A D^* F_n^* = (F_n D A^*)^*
        return self.FnD(A.conj().T).conj().T

    def spec(self):
        G, Hc = self.G, np.conj(self.H)
        lam, mu = self.lam, self.mu
        ent = lambda I, J: (G[I] @ Hc[J].T) / (lam[I][:, None] - mu[J][None, :])
        far_r = lambda I, zq: np.hstack([G[I, rr][:, None] / (lam[I][:, None] - zq[None, :])
                                         for rr in range(2)])
        far_c = lambda zq, J: np.vstack([Hc[J, rr][None, :] / (zq[:, None] - mu[J][None, :])
                                         for rr in range(2)])
        return KernelSpec(self.m, self.n, ent, lam, mu, far_r, far_c, cyclic=True,
                          dtype=complex, name='toeplitz_cauchy')


def make_conv_toeplitz(g, n, mode='valid'):
    """Toeplitz (tc, tr) of a 'valid' 1-D convolution with filter g (length 2w+1)."""
    w = (len(g) - 1) // 2
    m = n - 2 * w
    tc = np.zeros(m)
    tr = np.zeros(n)
    # row i: y_i = sum_j g(i + w - j) x_j, j in [i, i+2w]
    # T_{ij} = g[(i + w - j) + w] for |i+w-j| <= w -> depends on i-j  -> Toeplitz
    for d in range(-(n - 1), m):
        v = 0.0
        dd = d + w
        if -w <= dd <= w:
            v = g[dd + w]
        if d >= 0:
            tc[d] = v
        if d <= 0:
            tr[-d] = v
    return tc, tr
