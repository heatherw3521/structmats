function S = mp_exp_apps2(opts)
%MP_EXP_APPS2  Memo 2 applications (inner solves through H\b, minnorm, tikhonov).
%   BP   sparse-spike deconvolution (Ricker, stride 2) by basis pursuit:
%        IRLS (inner loop = diagonal-weight min-norm) and Douglas-Rachford
%        (inner step = min-norm projection H\b; H keeps its factors, so only
%        the first projection factors)
%   DB   Tikhonov deblurring, 1% noise: lambda sweep, discrepancy principle
%   MRI  non-Cartesian 1-D MRI: Cartesian k-space (diagonal) prior weights
if nargin < 1, opts = struct(); end
q = isfield(opts,'quick') && opts.quick;
rng(41);
n = 4096; if q, n = 1024; end
% ---------------- BP
K = mp_kernel_conv(n, 2, 'ricker', 2.5);
[L, rb, cb] = mp_tree(K.m, K.n, 64);
G = mp_hss_kernel(K, L, rb, cb, 1e-14); H = mp_hss_from_generators(G);
s = round(60*n/4096); pos = [];
while numel(pos) < s
    p = randi([3*K.w, n-3*K.w]);
    if all(abs(p - pos) >= 24), pos(end+1) = p; end %#ok<AGROW>
end
x0 = zeros(n,1); x0(pos) = sign(randn(s,1)) .* (0.5 + rand(s,1));
b = mp_conv_apply(K, x0);
xmn = H \ b;
x = xmn; ep = 1; hi = [];
t0 = tic;
for it = 1:60
    w = (x.^2 + ep^2).^(-0.25);
    x = mp_weighted_minnorm(H, b, w);
    rs = sort(abs(x), 'descend'); ep = max(min(ep, rs(s+1)/n), 1e-10);
    hi = [hi; it, toc(t0), mp_rel(x, x0), mp_rel(H*x, b)]; %#ok<AGROW>
    if mp_rel(x, x0) < 1e-10, break, end
end
S.BP.irls = hi; S.BP.x_irls = x;
gam = 0.05; z = xmn; hd = []; maxit = 400; if q, maxit = 100; end
t0 = tic;
proj = @(v) v + H \ (b - H*v);            % H keeps its factors: only the first call factors
soft = @(v, g) sign(v).*max(abs(v) - g, 0);
for it = 1:maxit
    xp = proj(z);
    z = z + soft(2*xp - z, gam) - xp;
    if mod(it, 10) == 0, hd = [hd; it, toc(t0), mp_rel(xp, x0), mp_rel(H*xp, b)]; end %#ok<AGROW>
end
S.BP.dr = hd; S.BP.x_dr = xp; S.BP.x0 = x0; S.BP.xmn = xmn; S.BP.m = K.m; S.BP.n = n; S.BP.s = s;
fprintf('BP m=%d n=%d s=%d | min-norm err %.2f | IRLS %d its err %.1e (%.1fs) | DR %d its err %.1e (%.1fs)\n', ...
    K.m, n, s, mp_rel(xmn, x0), size(hi,1), hi(end,3), hi(end,2), maxit, hd(end,3), hd(end,2));
% ---------------- DB
K = mp_kernel_conv(n, 2, 'gauss', 3);
[L, rb, cb] = mp_tree(K.m, K.n, 64);
H = mp_hss_from_generators(mp_hss_kernel(K, L, rb, cb, 1e-14));
t = linspace(0, 1, n)';
x0 = exp(-((t-0.2)/0.03).^2) + 0.6*(abs(t-0.5) < 0.08) + 0.8*exp(-((t-0.78)/0.01).^2) ...
     + 0.3*sin(6*pi*t).*(t > 0.85);
b0 = mp_conv_apply(K, x0); e = 1e-2*norm(b0)/sqrt(K.m)*randn(K.m,1); b = b0 + e;
lams = logspace(-5, 0.5, 23); rows = zeros(numel(lams), 4);
t0 = tic;
for j = 1:numel(lams)
    x = mp_tikhonov(H, b, lams(j));
    rows(j,:) = [lams(j), mp_rel(x, x0), norm(H*x - b), norm(x)];
end
tl = toc(t0)/numel(lams);
ok = rows(:,3) <= 1.01*norm(e);
lamdp = max(rows(ok,1)); [~, jb] = min(rows(:,2));
S.DB = struct('rows', rows, 'noise_norm', norm(e), 'lam_dp', lamdp, 'lam_best', rows(jb,1), ...
    'err_best', rows(jb,2), 'err_dp', mp_rel(mp_tikhonov(H, b, lamdp), x0), 'err_minnorm', mp_rel(H\b, x0), ...
    't_per_lambda', tl, 'x0', x0, 'x_dp', mp_tikhonov(H, b, lamdp));
fprintf('DB 1%% noise | min-norm err %.2e | best lam %.1e err %.3f | discrepancy lam %.1e err %.3f | %.2fs per lambda\n', ...
    S.DB.err_minnorm, S.DB.lam_best, S.DB.err_best, lamdp, S.DB.err_dp, tl);
% ---------------- MRI
N = 2048; if q, N = 512; end
tg = (-N/2:N/2-1)';
x0 = 1.0*(abs(tg) < 600*N/2048) + 0.5*(abs(tg - 150*N/2048) < 120*N/2048) ...
     - 0.4*(abs(tg + 260*N/2048) < 60*N/2048) + 0.3*exp(-((tg - 420*N/2048)/(25*N/2048)).^2);
m = round(0.45*N); u = sort(rand(m,1)); a = 0.12;
tot = a*log1p(0.5/a); qq = (u - 0.5)*2*tot;
xi = sign(qq).*a.*expm1(abs(qq)/a); xi = min(max(xi, -0.5+1e-9), 0.5-1e-9);
xp = sort(mod(-xi, 1));
K = mp_kernel_nudft(xp, N);
Av = exp(2i*pi*xp*tg.'); d = Av*x0;
Fu = exp(2i*pi*(0:N-1)'*tg.'/N)/sqrt(N);
[La, rba, cba] = mp_tree_aligned(xp, (0:N-1)'/N, 32);
H = mp_hss_from_generators(mp_hss_kernel(K, La, rba, cba, 1e-12));
l = (0:N-1)'; kf = l; kf(l >= N/2) = l(l >= N/2) - N;
wts = 1 + abs(kf)/8;
ymn = H \ d; yw = mp_weighted_minnorm(H, d, wts);
sg = 2e-2*norm(d)/sqrt(m); dn = d + sg*(randn(m,1) + 1i*randn(m,1))/sqrt(2);
lamg = logspace(-4, 0, 13); et = zeros(size(lamg));
for j = 1:numel(lamg), et(j) = mp_rel(Fu'*mp_tikhonov(H, dn, lamg(j), wts), x0); end
[eb, jb] = min(et);
S.MRI = struct('m', m, 'N', N, 'err_minnorm', mp_rel(Fu'*ymn, x0), 'err_weighted', mp_rel(Fu'*yw, x0), ...
    'err_minnorm_noisy', mp_rel(Fu'*(H\dn), x0), 'err_weighted_tik', eb, 'lam_best', lamg(jb), ...
    'x0', x0, 'x_mn', real(Fu'*ymn), 'x_w', real(Fu'*yw));
fprintf('MRI m=%d N=%d | image err: min-norm %.3f, k-space-weighted %.3f | noisy: min-norm %.3f, weighted Tikhonov %.3f\n', ...
    m, N, S.MRI.err_minnorm, S.MRI.err_weighted, S.MRI.err_minnorm_noisy, eb);
end
