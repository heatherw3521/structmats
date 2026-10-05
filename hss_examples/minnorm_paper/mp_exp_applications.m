function S = mp_exp_applications(opts)
%MP_EXP_APPLICATIONS  Application-driven tests of the min-norm solve
%   A1  band-limited (minimum-energy) interpolation from nonuniform samples
%       (NUDFT / seismic trace regularisation): jittered sampling and a gap;
%       index-based tree (hss_constructor rule) vs geometry-aligned tree
%   A2  1-D deconvolution, decimated Gaussian blur: exact and noisy data
%   A3  general Toeplitz through the Cauchy-like transform: HSS ranks of T vs C
%   A4  2-D blur (higher dimensions): HSS rank growth, lexicographic vs Morton
if nargin < 1, opts = struct(); end
q = isfield(opts,'quick') && opts.quick;
rng(11);
% ---------------- A1
N = 1024; if q, N = 256; end
k = (-N/2:N/2-1)';
c = (randn(N,1) + 1i*randn(N,1)) ./ (1 + abs(k)/40).^1.5;
c(N/2+2:end) = conj(c(N/2:-1:2)); c(N/2+1) = real(c(N/2+1)); c(1) = 0;   % real signal
f = @(x) real(exp(2i*pi*x(:)*k.') * c);
for cs = {'jitter', 'gap'}
    m = round(0.7*N);
    if strcmp(cs{1}, 'jitter')
        xs = sort(((0:m-1)' + 0.5 + 0.5*(rand(m,1)-0.5))/m);
    else
        m2 = round(1.15*m); xs = sort(((0:m2-1)' + 0.5 + 0.5*(rand(m2,1)-0.5))/m2);
        xs = xs(xs < 0.42 | xs > 0.52); xs = xs(1:min(m, numel(xs))); m = numel(xs);
    end
    K = mp_kernel_nudft(xs, N);
    fx = complex(f(xs));
    [L, rb, cb] = mp_tree(m, N, 32);
    [Gi, ii] = mp_hss_kernel(K, L, rb, cb, 1e-12);
    [La, rba, cba] = mp_tree_aligned(xs, (0:N-1)'/N, 32);
    [Ga, ia] = mp_hss_kernel(K, La, rba, cba, 1e-12);
    Hi = mp_hss_from_generators(Gi); Ha = mp_hss_from_generators(Ga);
    C = K.ent(1:m, 1:N); yr = mp_minnorm_dense(C, fx);
    yi = Hi \ fx; ya = Ha \ fx;
    Fu = exp(2i*pi*(0:N-1)'*k.'/N)/sqrt(N);
    cmn = Fu' * yi;
    tt = linspace(0, 1, 4000)';
    S.A1.(cs{1}) = struct('m', m, 'N', N, 'rank_index_tree', ii.maxrank, 'rank_aligned_tree', ia.maxrank, ...
        'err_index', mp_rel(yi, yr), 'err_aligned', mp_rel(ya, yr), 'norm_c_true', norm(c), ...
        'norm_c_minnorm', norm(cmn), 'xs', xs, 'fx', real(fx), 'tt', tt, 'ftrue', f(tt), ...
        'fmn', real(exp(2i*pi*tt*k.')*cmn), 'slack_ok_index', mp_slack_ok(Gi));
    fprintf('A1 %-6s m=%4d N=%4d  rank: index tree %d, aligned tree %d | err %.1e / %.1e | ||c_mn|| %.3f <= ||c|| %.3f\n', ...
        cs{1}, m, N, ii.maxrank, ia.maxrank, mp_rel(yi,yr), mp_rel(ya,yr), norm(cmn), norm(c));
end
% ---------------- A2
n = 2048; if q, n = 512; end
K = mp_kernel_conv(n, 2, 'gauss', 3);
[L, rb, cb] = mp_tree(K.m, K.n, 64);
H = mp_hss_from_generators(mp_hss_kernel(K, L, rb, cb, 1e-14));
t = linspace(0, 1, n)';
x0 = exp(-((t-0.2)/0.03).^2) + 0.6*(abs(t-0.5) < 0.08) + 0.8*exp(-((t-0.78)/0.01).^2) ...
     + 0.3*sin(6*pi*t).*(t > 0.85);
b = mp_conv_apply(K, x0);
A = K.ent(1:K.m, 1:K.n); s = svd(A);
xmn = H \ b; xr = mp_minnorm_dense(A, b);
e = 1e-3*norm(b)/sqrt(K.m)*randn(K.m,1);
xno = H \ (b + e);
S.A2 = struct('m', K.m, 'n', n, 'kappa', s(1)/s(end), 'err_vs_dense', mp_rel(xmn, xr), ...
    'dist_truth_exact', mp_rel(xmn, x0), 'dist_truth_noisy', mp_rel(xno, x0), 't', t, 'x0', x0, ...
    'xmn', xmn, 'xnoisy', xno, 'sv', s);
fprintf('A2 deconvolution %dx%d kappa=%.1e  err vs dense %.1e | ||x_mn - x0||/||x0||: exact %.2f, 0.1%% noise %.2e\n', ...
    K.m, n, s(1)/s(end), S.A2.err_vs_dense, S.A2.dist_truth_exact, S.A2.dist_truth_noisy);
% ---------------- A3
S.A3 = [];
for e = 9:ternary_(q, 10, 13)
    n = 2^e; m = n/2;
    tv = randn(m+n-1,1)/sqrt(n); tc = tv(m:-1:1); tr = tv(m:end); tc(1) = tr(1);
    K = mp_kernel_toeplitz(tc, tr);
    [L, rb, cb] = mp_tree(m, n, 64);
    [G, info] = mp_hss_kernel(K, L, rb, cb, 1e-12);
    H = mp_hss_from_generators(G);
    T = K.T();
    [GT, infoT] = mp_hss_kernel(struct('m', m, 'n', n, 'ent', @(I,J) T(I,J), 'rpos', complex((1:m)'), ...
        'cpos', complex((1:n)'), 'farr', [], 'farc', [], 'cyclic', false, 'isreal', true), L, rb, cb, 1e-12, ...
        struct('mode', 'full'));
    b = randn(m,1);
    t0 = tic; x = K.back(H \ K.fwd(b)); ts = toc(t0);
    S.A3 = [S.A3, struct('m', m, 'n', n, 'rank_C', info.maxrank, 'rank_T', infoT.maxrank, ...
        'err', mp_rel(x, mp_minnorm_dense(T, b)), 'imag_frac', norm(imag(x))/norm(x), 't_solve', ts)];
    fprintf('A3 Toeplitz %5dx%-5d  HSS rank of T %3d, of C %3d | err vs dense %.1e  (%.2fs)\n', ...
        m, n, infoT.maxrank, info.maxrank, S.A3(end).err, ts);
end
% ---------------- A4
S.A4 = [];
sig = 1.5; w = 5;
for n1 = ternary_(q, [16 32], [16 32 64])
    n2 = n1/2;
    [ii, jj] = ndgrid(0:n1-1, 0:n1-1); [oi, oj] = ndgrid(2*(0:n2-1)+0.5, 2*(0:n2-1)+0.5);
    for ord = {'lex', 'morton'}
        if strcmp(ord{1}, 'lex')
            pc = 1:n1^2; pr = 1:n2^2;
        else
            [~, pc] = sort(morton(ii(:), jj(:))); [~, pr] = sort(morton((oi(:)-0.5)/2, (oj(:)-0.5)/2));
        end
        ci = ii(pc); cj = jj(pc); ri = oi(pr); rj = oj(pr);
        ci = ci(:); cj = cj(:); ri = ri(:); rj = rj(:);
        ent = @(I,J) blur2(ri(I), rj(I), ci(J), cj(J), sig, w);
        K2 = struct('m', n2^2, 'n', n1^2, 'ent', ent, 'rpos', complex(zeros(n2^2,1)), ...
            'cpos', complex(zeros(n1^2,1)), 'farr', [], 'farc', [], 'cyclic', false, 'isreal', true);
        [L, rb, cb] = mp_tree(K2.m, K2.n, 32);
        [G, info] = mp_hss_kernel(K2, L, rb, cb, 1e-10, struct('mode', 'full'));
        H = mp_hss_from_generators(G);
        [b, xs] = mp_manufactured(H);
        [t, ~, x] = mp_timeit(@() H \ b, 1, 0);
        S.A4 = [S.A4, struct('n1', n1, 'N', n1^2, 'order', ord{1}, 'maxrank', info.maxrank, ...
            't', t, 'err', mp_rel(x, xs))];
        fprintf('A4 2-D blur %3dx%-3d (%s)  N=%5d  max HSS rank %3d  H\\b %.2fs  err %.1e\n', ...
            n1, n1, ord{1}, n1^2, info.maxrank, t, mp_rel(x, xs));
    end
end
end

function A = blur2(ri, rj, ci, cj, sig, w)
d2 = (ri - ci.').^2 + (rj - cj.').^2;
A = exp(-d2/(2*sig^2)); A(d2 > w^2) = 0;
end

function v = morton(i, j)
i = round(i(:)); j = round(j(:)); v = zeros(size(i)); b = 0;
while any(i > 0 | j > 0)
    v = v + bitand(j, 1)*2^(2*b) + bitand(i, 1)*2^(2*b+1);
    i = floor(i/2); j = floor(j/2); b = b + 1;
end
end

function v = ternary_(c, a, b)
if c, v = a; else, v = b; end
end
