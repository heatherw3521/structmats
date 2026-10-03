%MP_SELFTEST  Unit tests for the minnorm_paper helpers (small sizes, dense checks).
%   Run from hss_examples/minnorm_paper/ (adds lib/ and the repo root to the path).
%   Every line prints PASS/FAIL with the measured error.
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
addpath(fullfile(here, 'lib'));
if ~exist('OCTAVE_VERSION', 'builtin')
    addpath(fileparts(fileparts(here)));          % repo root (for @hss)
end
nfail = 0;
PF = {'PASS', 'FAIL'};
chk = @(name, err, tol) fprintf('%-58s %-4s  (%.2e)\n', name, PF{1 + (err > tol)}, err);

% ---- tree identical to the constructor's rule
[L, rb, cb] = mp_tree(320, 480, 15);
ok = (L == 4) && rb{5}(2) == 20 && cb{5}(2) == 30;
fprintf('%-58s %-4s\n', 'mp_tree(320,480,15): depth 4, first leaf 20x30', PF{1 + ~ok});

% ---- F1 generators -> hss object: full(H) == X*Y + blkdiag
G = mp_gen_lrbd(160, 240, 3, 15, 1);
H = mp_hss_from_generators(G);
A = full(H);
E = A - G.X*G.Y;
Lr = G.rb{G.L+1}; Lc = G.cb{G.L+1};
for i = 1:2^G.L, E(Lr(i)+1:Lr(i+1), Lc(i)+1:Lc(i+1)) = 0; end
chk('F1: off-leaf-diagonal part of full(H) equals X*Y', norm(E,'fro')/norm(A,'fro'), 1e-14);
x = randn(240,1);
chk('F1: H*x vs full(H)*x', norm(H*x - A*x)/norm(A*x), 1e-13);
chk('F1: H''*y vs full(H)''*y', norm(H'*x(1:160) - A'*x(1:160))/norm(A'*x(1:160)), 1e-13);
G2 = mp_hss_to_generators(H); H2 = mp_hss_from_generators(G2);
chk('generators round trip', norm(full(H2) - A,'fro')/norm(A,'fro'), 0);
b = randn(160,1);
xs = H \ b; xr = mp_minnorm_dense(A, b);
chk('F1: H\b vs dense min-norm', norm(xs - xr)/norm(xr), 1e-11);

% ---- F2 complex
G = mp_gen_random(128, 256, 6, 16, 2, struct('complex', true));
H = mp_hss_from_generators(G); A = full(H);
b = randn(128,1) + 1i*randn(128,1);
xs = H \ b; xr = mp_minnorm_dense(A, b);
chk('F2 complex: H\b vs dense min-norm', norm(xs - xr)/norm(xr), 1e-12);

% ---- kernel builders (proxy and full) vs dense kernel
K = mp_kernel_cauchy(256, 3);
[L, rb, cb] = mp_tree(K.m, K.n, 32);
Ad = K.ent(1:K.m, 1:K.n);
Gp = mp_hss_kernel(K, L, rb, cb, 1e-10, struct('mode','proxy'));
Gf = mp_hss_kernel(K, L, rb, cb, 1e-10, struct('mode','full'));
chk('F3 Cauchy, proxy build: ||H - A||_F/||A||_F', norm(full(mp_hss_from_generators(Gp)) - Ad,'fro')/norm(Ad,'fro'), 1e-8);
chk('F3 Cauchy, full build:  ||H - A||_F/||A||_F', norm(full(mp_hss_from_generators(Gf)) - Ad,'fro')/norm(Ad,'fro'), 1e-8);
K = mp_kernel_conv(600, 2, 'ricker', 2.5);
[L, rb, cb] = mp_tree(K.m, K.n, 32);
Ad = K.ent(1:K.m, 1:K.n);
Hc = mp_hss_from_generators(mp_hss_kernel(K, L, rb, cb, 1e-14));
xx = randn(K.n,1);
chk('F4 conv: proxy build vs dense', norm(full(Hc) - Ad,'fro')/norm(Ad,'fro'), 1e-13);
chk('F4 conv: mp_conv_apply vs dense', norm(mp_conv_apply(K, xx) - Ad*xx)/norm(Ad*xx), 1e-14);
rng(3); xs5 = sort(((0:299)' + 0.5 + 0.5*(rand(300,1)-0.5))/300);
K = mp_kernel_nudft(xs5, 384);
[L, rb, cb] = mp_tree(K.m, K.n, 24);
Ad = K.ent(1:K.m, 1:K.n);
kk = -192:191; Av = exp(2i*pi*xs5*kk); Fu = exp(2i*pi*(0:383)'*kk/384)/sqrt(384);
chk('F5 NUDFT: C = A F^*', norm(Ad - Av*Fu','fro')/norm(Ad,'fro'), 1e-12);
H5 = mp_hss_from_generators(mp_hss_kernel(K, L, rb, cb, 1e-12));
chk('F5 NUDFT: proxy build vs dense', norm(full(H5) - Ad,'fro')/norm(Ad,'fro'), 1e-10);
b = randn(300,1);
chk('F5 NUDFT: H\b vs dense min-norm', norm(H5\b - mp_minnorm_dense(full(H5), b))/norm(mp_minnorm_dense(full(H5), b)), 1e-11);
t = randn(128+256-1,1)/16; tc = t(128:-1:1); tr = t(128:end); tc(1) = tr(1);
K = mp_kernel_toeplitz(tc, tr);
T = K.T(); Ct = K.ent(1:128, 1:256);
Cfft = (fft(diag(exp(1i*pi*(0:255)'/256)) * (fft(T)/sqrt(128))')/sqrt(256))';
chk('F6 Toeplitz: Cauchy-like entries vs FFT transform', norm(Ct - Cfft,'fro')/norm(Cfft,'fro'), 1e-12);
[L, rb, cb] = mp_tree(128, 256, 16);
H6 = mp_hss_from_generators(mp_hss_kernel(K, L, rb, cb, 1e-12));
b = randn(128,1);
x6 = K.back(H6 \ K.fwd(b));
chk('F6 Toeplitz: back(H\fwd(b)) vs dense min-norm of T', norm(x6 - mp_minnorm_dense(T, b))/norm(mp_minnorm_dense(T,b)), 1e-9);

% ---- weighted min-norm and Tikhonov vs dense formulas
G = mp_gen_random(200, 400, 6, 20, 4); H = mp_hss_from_generators(G); A = full(H);
b = randn(200,1); w = 10.^(2*(rand(400,1)-0.5));
x = mp_weighted_minnorm(H, b, w);
xr = mp_minnorm_dense(A ./ w.', b) ./ w;
chk('weighted (diag): vs dense', norm(x - xr)/norm(xr), 1e-11);
cbL = G.cb{G.L+1}; Lb = cell(2^G.L,1); Lf = zeros(400);
for i = 1:2^G.L
    p = cbL(i+1)-cbL(i); M = randn(p); Lb{i} = chol(M*M' + p*eye(p));
    Lf(cbL(i)+1:cbL(i+1), cbL(i)+1:cbL(i+1)) = Lb{i};
end
x = mp_weighted_minnorm(H, b, [], Lb);
xr = Lf \ mp_minnorm_dense(A / Lf, b);
chk('weighted (leaf blocks): vs dense', norm(x - xr)/norm(xr), 1e-11);
for shape = {[200 400], [300 300], [400 200]}
    sz = shape{1};
    G = mp_gen_random(sz(1), sz(2), 6, 20, 5); H = mp_hss_from_generators(G); A = full(H);
    b = randn(sz(1),1);
    for lam = [1e-6, 1e-2, 1]
        [U_, S_, V_] = svd(A, 'econ'); s_ = diag(S_);
        xr = V_ * ((s_./(s_.^2 + lam^2)) .* (U_'*b));
        [x, r] = mp_tikhonov(H, b, lam);
        chk(sprintf('Tikhonov %dx%d lam=%.0e: vs SVD filter', sz(1), sz(2), lam), norm(x - xr)/norm(xr), 1e-10);
    end
end
chk('Tikhonov residual output r = Hx - b', norm(r - (A*x - b))/norm(A*x-b), 1e-10);
G = mp_gen_random(200, 300, 6, 20, 6); H = mp_hss_from_generators(G); A = full(H);
b = randn(200,1); w = 10.^(rand(300,1)-0.5); s = 10.^(rand(200,1)-0.5); lam = 0.1;
x = mp_tikhonov(H, b, lam, w, s);
xr = [s.*A; lam*diag(w)] \ [s.*b; zeros(300,1)];
chk('Tikhonov general form (S, W diagonal) vs stacked LS', norm(x - xr)/norm(xr), 1e-10);

% ---- proposed variant (proposed/mp_ulvminnormsolve_relaxed.m)
addpath(fullfile(here, 'proposed'));
G = mp_gen_random(256, 384, 12, 16, 7); H = mp_hss_from_generators(G); A = full(H);
b = randn(256,1); xr = mp_minnorm_dense(A, b);
fprintf('%-58s %-4s\n', 'proposed: test matrix violates Assumption 1', PF{1 + mp_slack_ok(G)});
chk('proposed (no Assumption 1): vs dense min-norm', norm(mp_ulvminnormsolve_relaxed(H, b) - xr)/norm(xr), 1e-11);
G = mp_gen_random(256, 256, 6, 16, 8); H = mp_hss_from_generators(G); A = full(H);
b = randn(256,1);
chk('proposed, square H: vs A\b', norm(mp_ulvminnormsolve_relaxed(H, b) - A\b)/norm(A\b), 1e-10);
K = mp_kernel_conv(300, 1, 'gauss', 2.0, 10);
[L, rb, cb] = mp_tree(K.m, K.n, 16);
H = mp_hss_from_generators(mp_hss_kernel(K, L, rb, cb, 1e-15)); A = full(H);
b = randn(K.m,1); bw = @(x) norm(b - A*x)/(norm(A)*norm(x) + norm(b));
chk('proposed, Gaussian blur (kappa ~ 1e8): backward error', bw(mp_ulvminnormsolve_relaxed(H, b)), 1e-15);
fprintf('%-58s INFO  (%.2e)\n', '  H\b, same matrix: backward error', bw(H\b));
