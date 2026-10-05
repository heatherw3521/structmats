function S = mp_exp_weighted_tikhonov(opts)
%MP_EXP_WEIGHTED_TIKHONOV  Weighted min-norm and Tikhonov accuracy/timing tests (all through H\b).
%   W1 weighted min-norm, L = diag(w), w_j = 10^(alpha*u_j), u ~ U[-1/2,1/2],
%      and leaf-block L_tau = chol(I + beta*K_tau) (Gaussian smoothing K_tau);
%      vs dense  x = L^{-1} (H L^{-1})^+ b
%   W2 large-scale weighted, manufactured  x* = L^{-1} (H L^{-1})^* z
%   T1 Tikhonov lambda sweep for wide / square / tall H vs the SVD filter
%      solution; general form min ||S(Hx-b)||^2 + lam^2 ||W x||^2
%   T2 large-scale Tikhonov, manufactured  x* = H^* z, b = H x* + lam^2 z,
%      and time vs plain H\b
if nargin < 1, opts = struct(); end
q = isfield(opts,'quick') && opts.quick;
maxexp = 18; if isfield(opts,'maxexp'), maxexp = opts.maxexp; end
if q, maxexp = 12; end
rng(31);
% ---------------- W1
S.W1 = [];
m = 1024; n = 2048; if q, m = 256; n = 512; end
for fam = {'F2', 'F3', 'F4'}
    switch fam{1}
        case 'F2', G = mp_gen_random(m, n, 10, 32);
        case 'F3', K = mp_kernel_cauchy(round(n/3)+1, 3); [L,rb,cb] = mp_tree(K.m,K.n,32); G = mp_hss_kernel(K,L,rb,cb,1e-12);
        case 'F4', K = mp_kernel_conv(n, 2, 'gauss', 2); [L,rb,cb] = mp_tree(K.m,K.n,32); G = mp_hss_kernel(K,L,rb,cb,1e-14);
    end
    H = mp_hss_from_generators(G); A = full(H); [mm, nn] = size(A);
    b = randn(mm,1); u = rand(nn,1) - 0.5;
    for alpha = [0 2 4 6 8 10 12]
        w = 10.^(alpha*u);
        x = mp_weighted_minnorm(H, b, w);
        Aw = A ./ w.'; xr = mp_minnorm_dense(Aw, b) ./ w; s = svd(Aw);
        S.W1 = [S.W1, struct('family', fam{1}, 'kind', 'diag', 'param', alpha, 'kappa', s(1)/s(end), ...
            'err', mp_rel(x, xr), 'res', mp_rel(A*x, b))];
        fprintf('W1 %s diag alpha=%2d kappa(HL^-1)=%.1e err %.1e res %.1e\n', fam{1}, alpha, s(1)/s(end), S.W1(end).err, S.W1(end).res);
    end
    cbL = G.cb{G.L+1};
    for beta = [1 1e2 1e4]
        Lb = cell(2^G.L,1); Lf = zeros(nn);
        for i = 1:2^G.L
            p = cbL(i+1)-cbL(i); tt = (0:p-1)';
            Kt = exp(-(tt - tt.').^2/(2*3^2));
            Lb{i} = chol(eye(p) + beta*Kt);
            Lf(cbL(i)+1:cbL(i+1), cbL(i)+1:cbL(i+1)) = Lb{i};
        end
        x = mp_weighted_minnorm(H, b, [], Lb);
        ALi = A / Lf; xr = Lf \ mp_minnorm_dense(ALi, b); s = svd(ALi);
        S.W1 = [S.W1, struct('family', fam{1}, 'kind', 'block', 'param', beta, 'kappa', s(1)/s(end), ...
            'err', mp_rel(x, xr), 'res', mp_rel(A*x, b))];
        fprintf('W1 %s block beta=%g kappa=%.1e err %.1e\n', fam{1}, beta, s(1)/s(end), S.W1(end).err);
    end
end
% ---------------- W2
S.W2 = [];
for e = 12:maxexp
    n = 2^e;
    G = mp_gen_random(n/2, n, 10, 64); H = mp_hss_from_generators(G);
    w = 10.^(6*(rand(n,1) - 0.5));
    Hw = mp_hss_from_generators(mp_gen_scale(G, [], 1./w));
    z = randn(n/2,1); xs = (Hw'*z) ./ w; b = H*xs;
    [t, ~, x] = mp_timeit(@() mp_weighted_minnorm(H, b, w), 1, 0);
    S.W2 = [S.W2, struct('n', n, 't', t, 'err', mp_rel(x, xs), 'res', mp_rel(H*x, b))];
    fprintf('W2 n=%8d  %.2fs  err %.1e  res %.1e\n', n, t, mp_rel(x, xs), mp_rel(H*x, b));
end
% ---------------- T1
S.T1 = [];
shapes = {[1024 2048], [1536 1536], [2048 1024]};
if q, shapes = {[256 512], [384 384], [512 256]}; end
for sh = shapes
    sz = sh{1};
    G = mp_gen_random(sz(1), sz(2), 10, 32); H = mp_hss_from_generators(G); A = full(H);
    [U_, S_, V_] = svd(A, 'econ'); s = diag(S_); b = randn(sz(1),1);
    for lr = -8:2
        lam = 10^lr * s(1);
        x = mp_tikhonov(H, b, lam);
        xr = V_ * ((s ./ (s.^2 + lam^2)) .* (U_'*b));
        S.T1 = [S.T1, struct('m', sz(1), 'n', sz(2), 'lam_rel', 10^lr, 'err', mp_rel(x, xr))];
        fprintf('T1 %dx%d lam/||A||=%.0e  err %.1e\n', sz(1), sz(2), 10^lr, mp_rel(x, xr));
    end
end
G = mp_gen_random(m, round(1.5*m), 10, 32); H = mp_hss_from_generators(G); A = full(H);
[mm, nn] = size(A); b = randn(mm,1);
for lam = [1e-4 1e-2 1]
    w = 10.^(2*(rand(nn,1)-0.5)); sr = 10.^(2*(rand(mm,1)-0.5));
    x = mp_tikhonov(H, b, lam, w, sr);
    xr = [sr.*A; lam*diag(w)] \ [sr.*b; zeros(nn,1)];
    S.T1 = [S.T1, struct('m', mm, 'n', nn, 'lam_rel', -lam, 'err', mp_rel(x, xr))];
    fprintf('T1 general form lam=%.0e  err %.1e\n', lam, mp_rel(x, xr));
end
% ---------------- T2
S.T2 = [];
for e = 10:maxexp
    n = 2^e; lam = 1e-2;
    for shape = {'wide', 'tall'}
        if strcmp(shape{1}, 'wide'), mm = n/2; else, mm = 2*n; end
        G = mp_gen_random(mm, n, 10, 64); H = mp_hss_from_generators(G);
        z = randn(mm,1); xs = H'*z; b = H*xs + lam^2*z;
        [t, ~, x] = mp_timeit(@() mp_tikhonov(H, b, lam), 1, 0);
        tm = NaN;
        if strcmp(shape{1}, 'wide'), tm = mp_timeit(@() H \ b, 1, 0); end
        S.T2 = [S.T2, struct('n', n, 'm', mm, 'shape', shape{1}, 't_tikhonov', t, 't_minnorm', tm, 'err', mp_rel(x, xs))];
        fprintf('T2 %s %7dx%-7d  Tikhonov %.2fs  (H\\b %.2fs)  err %.1e\n', shape{1}, mm, n, t, tm, mp_rel(x, xs));
    end
end
end
