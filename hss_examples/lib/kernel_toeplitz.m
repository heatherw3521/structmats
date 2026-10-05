function K = kernel_toeplitz(tc, tr, phi)
%KERNEL_TOEPLITZ  Family F6, Cauchy-like transform of a rectangular Toeplitz
%   matrix T (m x n), T(i,j) = t_{i-j}, tc = T(:,1), tr = T(1,:).
%   With Z_theta = cyclic down-shift carrying theta = e^{i phi} in its corner,
%     Z_1^{(m)} T - T Z_theta^{(n)} = e_1 r.' + c e_n.'          (rank <= 2)
%   and the unitary maps F_m (F(j,k) = e^{-2 pi i jk/m}/sqrt(m)), D = diag(theta^{k/n}),
%     C = F_m T D^* F_n^*,   Lam C - C Mu = G H^*,
%     lam_j = e^{-2 pi i j/m},  mu_l = theta^{1/n} e^{-2 pi i l/n},
%     G = F_m [e_1, c],   H = F_n D [conj(r), e_n].
%   T x = b  <=>  C y = F_m b  with  y = F_n D x  (norm preserving), so the
%   minimum-norm solutions correspond.
%   K.fwd(b) = F_m b,  K.back(y) = D^* F_n^* y.  phi = pi (repo default, theta = -1).
if nargin < 3, phi = pi; end
tc = tc(:); tr = tr(:);
m = numel(tc); n = numel(tr);
theta = exp(1i*phi);
Tent = @(i,j) tval(tc, tr, i - j);                % 0-based i, j
lastrow = arrayfun(@(jj) Tent(m-1, jj), (0:n-1)');
r = lastrow - [tr(2:end); theta*tr(1)];
lastcol = arrayfun(@(ii) Tent(ii, n-1), (0:m-1)');
c = zeros(m,1); c(2:end) = lastcol(1:end-1) - theta*tc(2:end);
d = exp(1i*phi*(0:n-1)'/n);
Fm = @(v) fft(v)/sqrt(m);
FnD = @(v) fft(d .* v)/sqrt(n);
e1 = zeros(m,1); e1(1) = 1; en = zeros(n,1); en(end) = 1;
Gg = [Fm(e1), Fm(c)];
Hh = [FnD(conj(r)), FnD(en)];
lam = exp(-2i*pi*(0:m-1)'/m);
mu = exp(1i*phi/n) * exp(-2i*pi*(0:n-1)'/n);
Hc = conj(Hh);
K.m = m; K.n = n;
K.ent = @(I,J) (Gg(I,:) * Hc(J,:).') ./ (lam(I) - mu(J).');
K.rpos = lam; K.cpos = mu;
K.farr = @(I,zq) [Gg(I,1) ./ (lam(I) - zq.'), Gg(I,2) ./ (lam(I) - zq.')];
K.farc = @(zq,J) [Hc(J,1).' ./ (zq - mu(J).'); Hc(J,2).' ./ (zq - mu(J).')];
K.cyclic = true; K.isreal = false;
K.fwd = Fm;
K.back = @(y) conj(d) .* ifft(y) * sqrt(n);
K.T = @() toeplitz(tc, tr);
K.name = 'toeplitz_cauchy';
end

function v = tval(tc, tr, dd)
if dd >= 0, v = tc(dd+1); else, v = tr(-dd+1); end
end
