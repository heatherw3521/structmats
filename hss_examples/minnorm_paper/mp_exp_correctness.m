function R = mp_exp_correctness(opts)
%MP_EXP_CORRECTNESS  E1 (memo 1, Section 7.3, Table 2): H\b against the dense minimum-norm
%   solution, for every test family at a size where dense references fit.
%   Records kappa(H), the residual, the error against the dense minimum-norm
%   solution of the HSS matrix itself (solver error alone), the error against
%   the minimum-norm solution of the original matrix A (adds the compression
%   error; for the exact families F1 and F2, A is the HSS matrix), the
%   null-space fraction of x, and whether the leaves satisfy the slack
%   condition of Assumption 1 (the solver does not need it).
%   Leaves are small (16-32 rows) so that every tree has 4-7 levels.
%   Writes results/memo1_correctness.txt.
if nargin < 1, opts = struct(); end
q = isfield(opts, 'quick') && opts.quick;
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
rng(2026);
R = struct('name', {}, 'family', {}, 'm', {}, 'n', {}, 'L', {}, 'maxrank', {}, 'kappa', {}, ...
           'res', {}, 'err_H', {}, 'err_A', {}, 'null', {}, 'slack_ok', {}, 'time', {});
add = @(R, r) [R, r];
sz = @(a, b) ternary_(q, b, a);
% ---------------- F1: low rank plus block diagonal (exact, generators)
R = add(R, onecase('F1', 'F1', mp_gen_lrbd(320, 480, 3, 16), []));
R = add(R, onecase('F1', 'F1', mp_gen_lrbd(sz(2048,512), sz(4096,1024), 3, 16), []));
% ---------------- F2: random HSS (exact, generators)
R = add(R, onecase('F2', 'F2', mp_gen_random(sz(1024,256), sz(2048,512), 10, 32), []));
R = add(R, onecase('F2 complex', 'F2', mp_gen_random(sz(1024,256), sz(2048,512), 10, 32, [], struct('complex', true)), []));
R = add(R, onecase('F2 m/n = 0.9', 'F2', mp_gen_random(sz(1800,450), sz(2000,500), 8, 32), []));
% ---------------- F3: interlaced Cauchy (nested IDs)
K = mp_kernel_cauchy(sz(1024,256), 3);
[L, rb, cb] = mp_tree(K.m, K.n, 32);
R = add(R, onecase('F3', 'F3', mp_hss_kernel(K, L, rb, cb, 1e-10), K.ent(1:K.m, 1:K.n)));
% ---------------- F4: decimated convolution (exact, banded)
cf = {'gauss', 2, 2; 'ricker', 3, 2; 'gauss', 3, 3};
for c = 1:size(cf,1)
    K = mp_kernel_conv(sz(4096,1024), cf{c,3}, cf{c,1}, cf{c,2});
    [L, rb, cb] = mp_tree(K.m, K.n, 32);
    R = add(R, onecase(sprintf('F4 %s sigma=%g s=%d', cf{c,1}, cf{c,2}, cf{c,3}), 'F4', ...
        mp_hss_kernel(K, L, rb, cb, 1e-13), K.ent(1:K.m, 1:K.n)));
end
% ---------------- F5: NUDFT Cauchy-like (nested IDs)
mm = sz(1536,384); NN = sz(2048,512);
xs = sort(((0:mm-1)' + 0.5 + 0.5*(rand(mm,1) - 0.5))/mm);
K = mp_kernel_nudft(xs, NN);
[L, rb, cb] = mp_tree(K.m, K.n, 32);
R = add(R, onecase('F5', 'F5', mp_hss_kernel(K, L, rb, cb, 1e-12), K.ent(1:K.m, 1:K.n)));
% ---------------- F6: random Toeplitz through the Cauchy-like transform
m6 = sz(1024,256); n6 = 2*m6;
t = randn(m6+n6-1,1)/sqrt(n6); tc = t(m6:-1:1); tr = t(m6:end); tc(1) = tr(1);
K = mp_kernel_toeplitz(tc, tr);
[L, rb, cb] = mp_tree(m6, n6, 32);
G6 = mp_hss_kernel(K, L, rb, cb, 1e-12);
r6 = onecase('F6', 'F6', G6, K.ent(1:m6, 1:n6));
H6 = mp_hss_from_generators(G6); T = K.T(); b = randn(m6,1);
x = K.back(H6 \ K.fwd(b));
r6.err_A = mp_rel(x, mp_minnorm_dense(T, b));     % end to end, in the Toeplitz variables
R = add(R, r6);
for k = 1:numel(R)
    r = R(k);
    fprintf('%-28s %5dx%-6d L=%d r=%3d kappa=%.1e | res %.1e err(H) %.1e err(A) %.1e null %.1e slack_ok=%d  %.3fs\n', ...
        r.name, r.m, r.n, r.L, r.maxrank, r.kappa, r.res, r.err_H, r.err_A, r.null, r.slack_ok, r.time);
end
tag = ''; if q, tag = '_quick'; end
mp_write_table(fullfile(here, 'results', ['memo1_correctness' tag '.txt']), R, ...
    {'name', 'family', 'm', 'n', 'L', 'maxrank', 'kappa', 'res', 'err_H', 'err_A', 'null', 'slack_ok', 'time'}, ...
    'E1 correctness (memo 1, Table 2): err_H vs dense min-norm of H, err_A vs dense min-norm of the original matrix');
end

function r = onecase(name, family, G, A)
H = mp_hss_from_generators(G);
Ah = full(H);
if isempty(A), A = Ah; end                % exact families: the HSS matrix is the matrix
cplx = ~isreal(Ah);
b = randn(size(Ah,1),1); if cplx, b = b + 1i*randn(size(b)); end
t0 = tic; x = H \ b; t = toc(t0);
xr = mp_minnorm_dense(Ah, b);
s = svd(Ah);
mr = 0;
for l = 2:G.L+1
    for i = 1:numel(G.U{l}), mr = max([mr, size(G.U{l}{i},2), size(G.V{l}{i},2)]); end
end
r.name = name; r.family = family; r.m = size(Ah,1); r.n = size(Ah,2); r.L = G.L; r.maxrank = mr;
r.kappa = s(1)/s(end); r.res = norm(Ah*x - b)/norm(b); r.err_H = mp_rel(x, xr);
r.err_A = mp_rel(x, mp_minnorm_dense(A, b));
r.null = mp_nullfrac(Ah, x); r.slack_ok = mp_slack_ok(G); r.time = t;
end

function v = ternary_(c, a, b)
if c, v = a; else, v = b; end
end
