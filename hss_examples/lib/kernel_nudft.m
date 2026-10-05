function K = kernel_nudft(xs, N)
%KERNEL_NUDFT  Family F5, Cauchy-like form of the type-II NUDFT.
%   A(j,k) = exp(2 pi i k x_j), k = -N/2..N/2-1, x_j in [0,1) sorted, N even.
%   F unitary DFT, (F c)_l = sum_k c_k exp(2 pi i k l/N)/sqrt(N).  Then
%   C = A F^* has C(j,l) = D_N(x_j - l/N)/sqrt(N),
%   D_N(u) = exp(-i pi u) sin(pi N u)/sin(pi u)   (Dirichlet kernel), and
%   C(j,l) = a_j b_l/(z_j - w_l), z = e^{2 pi i x}, w_l = e^{2 pi i l/N}.
%   min ||c|| s.t. A c = f  <=>  min ||y|| s.t. C y = f  (y = F c, unitary).
xs = xs(:); m = numel(xs);
grid = (0:N-1)'/N;
z = exp(2i*pi*xs); wl = exp(2i*pi*grid);
a = 2i*sin(pi*N*xs)/sqrt(N);
b = (-1).^(0:N-1)' .* exp(2i*pi*grid);
K.m = m; K.n = N; K.xs = xs;
K.ent = @(I,J) dirichlet(xs(I), grid(J), N);
K.rpos = z; K.cpos = wl;
K.farr = @(I,zq) a(I) ./ (z(I) - zq.');
K.farc = @(zq,J) b(J).' ./ (zq - wl(J).');
K.cyclic = true; K.isreal = false;
K.name = 'nudft_cauchy';
end

function C = dirichlet(x, g, N)
u = x - g.';
num = sin(pi*N*u); den = sin(pi*u);
C = exp(-1i*pi*u) .* num ./ (den*sqrt(N));
small = abs(den) < 1e-14;
C(small) = sqrt(N);
end
