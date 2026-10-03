function nfail = mp_test_solver()
%MP_TEST_SOLVER  Regression checks for the ULV solver behind H\b (wide and
%   square H, @hss/private/hss_ulvminnormsolve.m), the factors H keeps
%   between solves, and the methods minnorm (weighted) and tikhonov.  Prints
%   PASS/FAIL per check, like mp_selftest, and returns the number of
%   failures.  Run from hss_examples/minnorm_paper/.
%
%   Every solve is compared with a dense reference at moderate size:
%     err    ||x - x_ref|| / ||x_ref||, x_ref = dense minimum-norm solution
%            (QR of H'), allowed up to max(1e-11, 100*kappa*eps);
%     bwd    normwise backward error ||b - A x|| / (||A|| ||x|| + ||b||), <= 1e-14;
%     levels number of levels in the factorization; must equal the tree depth
%            (the legacy solver merged levels when l + m > n at a leaf).
%   Covered: families F1-F6; leaves that violate the legacy slack assumption;
%   a leaf level with t = 0 decoupled rows; a depth-6 tree; an uneven,
%   geometry-aligned tree; matrices from the hss() constructor (the cases of
%   hss_examples/tests/test_hss.m); ill-conditioned Gaussian blur; complex
%   data; square systems; several right-hand sides; stored-factor reuse and
%   invalidation; weighted minimum norm and Tikhonov with diagonal and
%   leaf-block weights; and the two named errors (rankDeficient, improperRanks).
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
addpath(fullfile(here, 'lib'));
if ~exist('OCTAVE_VERSION', 'builtin')
    addpath(fileparts(fileparts(here)));          % repo root (for @hss)
end
nfail = 0;
fprintf('%-44s %-11s %-9s %-9s %-7s\n', 'case', 'm x n', 'err', 'bwd', 'levels');

% ---------------------------------------------------------------- families
C = {};
C{end+1} = {'F1 LR + block diagonal, k=3', mp_gen_lrbd(512, 1024, 3, 16, 1)};
C{end+1} = {'F2 random HSS, k=10', mp_gen_random(512, 1024, 10, 32, 2)};
C{end+1} = {'F2 complex, k=10', mp_gen_random(512, 1024, 10, 32, 3, struct('complex', true))};
C{end+1} = {'F2 m/n=0.50, k=16', mp_gen_random(512, 1024, 16, 64, 4)};
C{end+1} = {'F2 m/n=0.90, k=16 (legacy slack fails)', mp_gen_random(922, 1024, 16, 64, 4)};
C{end+1} = {'F2 m/n=0.95, k=16 (legacy slack fails)', mp_gen_random(973, 1024, 16, 64, 4)};
C{end+1} = {'F2 16x24 leaves, k=12 (l+m > n)', mp_gen_random(256, 384, 12, 16, 7)};
C{end+1} = {'F2 k = leaf rows (t = 0 at leaves)', mp_gen_random(256, 512, 16, 16, 5)};
C{end+1} = {'F2 depth-6 tree, k=8', mp_gen_random(1024, 2048, 8, 16, 6)};
K = mp_kernel_cauchy(342, 3); [L, rb, cb] = mp_tree(K.m, K.n, 32);
C{end+1} = {'F3 interlaced Cauchy', mp_hss_kernel(K, L, rb, cb, 1e-12)};
K = mp_kernel_conv(2048, 2, 'ricker', 3); [L, rb, cb] = mp_tree(K.m, K.n, 32);
C{end+1} = {'F4 Ricker convolution, stride 2', mp_hss_kernel(K, L, rb, cb, 1e-14)};
rng(11); x5 = sort(((0:767)' + 0.5 + 0.5*(rand(768,1)-0.5))/768);
K = mp_kernel_nudft(x5, 1024); [L, rb, cb] = mp_tree(K.m, K.n, 32);
C{end+1} = {'F5 NUDFT Cauchy-like (legacy slack fails)', mp_hss_kernel(K, L, rb, cb, 1e-12)};
t6 = randn(512+1024-1, 1)/16; tc = t6(512:-1:1); tr = t6(512:end); tc(1) = tr(1);
K = mp_kernel_toeplitz(tc, tr); [L, rb, cb] = mp_tree(512, 1024, 32);
C{end+1} = {'F6 Toeplitz -> Cauchy-like', mp_hss_kernel(K, L, rb, cb, 1e-12)};
xg = sort(rand(700, 1)); xg = xg(xg < 0.42 | xg > 0.52); N = 1024;
K = mp_kernel_nudft(xg, N); [L, rb, cb] = mp_tree_aligned(mod(-xg, 1), (0:N-1)'/N, 32);
C{end+1} = {'F5 with a gap, geometry-aligned tree', mp_hss_kernel(K, L, rb, cb, 1e-12)};
for sig = [1.6 2.0]
    K = mp_kernel_conv(300, 1, 'gauss', sig, ceil(5*sig)); [L, rb, cb] = mp_tree(K.m, K.n, 16);
    C{end+1} = {sprintf('Gaussian blur sigma=%.1f (ill-conditioned)', sig), ...
                mp_hss_kernel(K, L, rb, cb, 1e-15)}; %#ok<AGROW>
end
for c = 1:numel(C)
    H = mp_hss_from_generators(C{c}{2});
    nfail = nfail + solvecase(C{c}{1}, H, full(H), 20 + c);
end

% ---------------------------------------- matrices from the hss() constructor
rng(1); A = nested_lowrank(32, 64, 16, 2);
nfail = nfail + solvecase('hss(): nestedLowRank 32x64, k=2', hss(A, 'blocksize', 16, 'k', 2), [], 41);
rng(1); A = nested_lowrank(16, 20, 8, 6);
nfail = nfail + solvecase('hss(): 16x20, rank 6 > leaf slack', hss(A, 'blocksize', 8, 'tol', 1e-10), [], 42);
rng(2); A = nested_lowrank(256, 300, 16, 6);
nfail = nfail + solvecase('hss(): nestedLowRank 256x300', hss(A, 'blocksize', 16, 'tol', 1e-12), [], 43);

% ---------------------------------------------------------- square, via H\b
try
    H = mp_hss_from_generators(mp_gen_random(1024, 1024, 10, 32, 9)); A = full(H);
    rng(44); b = randn(1024, 1);
    x = H \ b;
    e = mp_rel(x, A\b); lv = depth_of(H.factorcache.ulv);
    ok = e < 1e-10 && lv == H.levelcount;
    fprintf('%-44s %-11s %-9.1e %-9s %d/%-5d %s\n', 'square H\b (vs A\b)', '1024x1024', e, '', ...
            lv, H.levelcount, pf(ok));
    nfail = nfail + ~ok;
    n = 400; A = 1 ./ ((1:n)'/n + (1:n)/n + 0.5) + n*eye(n);
    H = hss(A, 'blocksize', 25, 'tol', 1e-12); rng(48); b = randn(n, 1);
    nfail = nfail + chk('square hss() Cauchy + n*I: H\b vs A\b', mp_rel(H \ b, A \ b), 1e-10);
catch err
    nfail = nfail + report_error('square H\b', err);
end

% --------------------------------------- several right-hand sides, factor reuse
try
    H = mp_hss_from_generators(mp_gen_random(512, 1024, 10, 32, 8)); A = full(H);
    rng(45); B = randn(512, 4);
    Xr = mp_minnorm_dense(A, B);
    X = H \ B;
    Xc = zeros(1024, 4); for j = 1:4, Xc(:, j) = H \ B(:, j); end
    nfail = nfail + chk('H\B, 4 columns: vs column by column', mp_rel(X, Xc), 1e-12);
    nfail = nfail + chk('H\B, 4 columns: vs dense min-norm', mp_rel(X, Xr), 1e-12);
    % the factors H keeps after its first solve
    nfail = nfail + chk('stored factors: present after H\b', double(isempty(H.factorcache.ulv)), 0);
    nfail = nfail + chk('stored factors: one entry per level', abs(depth_of(H.factorcache.ulv) - H.levelcount), 0);
    x1 = H \ B(:, 1);                    % from the stored factors
    clearfactors(H);
    nfail = nfail + chk('clearfactors empties the store', double(~isempty(H.factorcache.ulv)), 0);
    nfail = nfail + chk('stored-factor solve equals a fresh factorization', mp_rel(x1, H \ B(:, 1)), 0);
    H2 = H;                               % a copy shares the factors ...
    nfail = nfail + chk('copy H2 = H shares the stored factors', double(isempty(H2.factorcache.ulv)), 0);
    H2.A11.D = 2*H2.A11.D;                % ... until it is modified
    nfail = nfail + chk('modifying the copy drops its factors', double(~isempty(H2.factorcache.ulv)), 0);
    nfail = nfail + chk('modified copy: H2\b vs dense', mp_rel(H2 \ B(:, 1), mp_minnorm_dense(full(H2), B(:, 1))), 1e-12);
    nfail = nfail + chk('original keeps its own factors and answer', mp_rel(H \ B(:, 1), Xr(:, 1)), 1e-12);
    % projection onto {x : Hx = b}: P(v) = v + H^+(b - H v) is idempotent and feasible
    rng(46); v = randn(1024, 1); bb = B(:, 3);
    Pv = v + H \ (bb - H*v); PPv = Pv + H \ (bb - H*Pv);
    nfail = nfail + chk('projection with stored factors: P(P(v)) = P(v)', mp_rel(PPv, Pv), 1e-12);
    nfail = nfail + chk('projection with stored factors: H P(v) = b', mp_rel(H*Pv, bb), 1e-12);
    nfail = nfail + chk('projection: P(v) - v orthogonal to null(H)', ...
                        norm(null_part(A, Pv - v))/norm(Pv - v), 1e-10);
catch err
    nfail = nfail + report_error('several right-hand sides / factor reuse', err);
end

% --------------------------------- weighted minimum norm and Tikhonov (methods)
for cplx = [false true]
    try
        G = mp_gen_random(256, 512, 8, 32, 12, struct('complex', cplx)); H = mp_hss_from_generators(G);
        A = full(H); [m, n] = size(A); rng(49); b = randn(m, 2) + 1i*cplx*randn(m, 2);
        cb = G.cb{end}; rb = G.rb{end}; nl = numel(cb) - 1;
        w = 10.^(2*(rand(n, 1) - 0.5)); s = 1 + rand(m, 1);
        Lb = cell(nl, 1); Sb = cell(nl, 1);
        for i = 1:nl
            p = cb(i+1) - cb(i); q = rb(i+1) - rb(i);
            Lb{i} = 3*eye(p) + 0.5*randn(p); Sb{i} = eye(q) + 0.3*randn(q);
        end
        Lf = blkdiag(Lb{:}); Sf = blkdiag(Sb{:}); lam = 0.3; tag = sprintf(' (complex=%d)', cplx);
        xr = diag(1./w)*mp_minnorm_dense(A*diag(1./w), b);
        nfail = nfail + chk(['minnorm Weight=vector vs dense' tag], mp_rel(minnorm(H, b, 'Weight', w), xr), 1e-11);
        xr = Lf \ mp_minnorm_dense(A/Lf, b);
        x = minnorm(H, b, 'Weight', Lb);
        nfail = nfail + chk(['minnorm Weight=leaf blocks vs dense' tag], mp_rel(x, xr), 1e-11);
        nfail = nfail + chk(['minnorm repeat (stored factors) identical' tag], mp_rel(minnorm(H, b, 'Weight', Lb), x), 0);
        xr = [Sf*A; lam*Lf] \ [Sf*b; zeros(n, 2)];
        [x, r] = tikhonov(H, b, lam, 'Weight', Lb, 'DataWeight', Sb);
        nfail = nfail + chk(['tikhonov L, S leaf blocks vs dense' tag], mp_rel(x, xr), 1e-10);
        nfail = nfail + chk(['tikhonov residual r = S(Hx - b)' tag], mp_rel(r, Sf*(A*x - b)), 1e-10);
        xr = [diag(s)*A; lam*diag(w)] \ [diag(s)*b; zeros(n, 2)];
        nfail = nfail + chk(['tikhonov L, S vectors vs dense' tag], ...
                            mp_rel(tikhonov(H, b, lam, 'Weight', w, 'DataWeight', s), xr), 1e-10);
    catch err
        nfail = nfail + report_error(sprintf('weighted / Tikhonov (complex=%d)', cplx), err);
    end
end
try
    G = mp_gen_random(512, 256, 8, 32, 13); H = mp_hss_from_generators(G); A = full(H);
    rng(50); b = randn(512, 1); cb = G.cb{end}; nl = numel(cb) - 1; Lb = cell(nl, 1);
    for i = 1:nl, p = cb(i+1) - cb(i); Lb{i} = triu(randn(p)) + 4*eye(p); end
    lam = 1e-3; xr = [A; lam*blkdiag(Lb{:})] \ [b; zeros(256, 1)];
    nfail = nfail + chk('tikhonov, tall H, L leaf blocks vs dense', mp_rel(tikhonov(H, b, lam, 'Weight', Lb), xr), 1e-9);
catch err
    nfail = nfail + report_error('tikhonov, tall H', err);
end

% ---------------------------------------------------------------- named errors
G = struct('L', 1, 'rb', {{[0 20], [0 12 20]}}, 'cb', {{[0 24], [0 8 24]}});
G.D = {randn(12, 8), randn(8, 16)};
G.U = {{}, {orth(randn(12, 2)), orth(randn(8, 2))}};
G.V = {{}, {orth(randn(8, 2)), orth(randn(16, 2))}};
G.B12 = {{randn(2)}}; G.B21 = {{randn(2)}};
H = mp_hss_from_generators(G);           % leaf 1: 10 uncoupled rows on 8 columns
nfail = nfail + expect_error('rank-deficient H -> named error', @() H \ randn(20, 1), ...
                             'hss_ulvminnormsolve:rankDeficient');
G = struct('L', 1, 'rb', {{[0 8], [0 4 8]}}, 'cb', {{[0 24], [0 12 24]}});
G.D = {randn(4, 12), randn(4, 12)};
G.U = {{}, {randn(4, 6), randn(4, 6)}};  % row basis of rank 6 on a 4-row leaf
G.V = {{}, {orth(randn(12, 6)), orth(randn(12, 6))}};
G.B12 = {{randn(6)}}; G.B21 = {{randn(6)}};
H = mp_hss_from_generators(G);
nfail = nfail + expect_error('row basis wider than leaf -> named error', @() H \ randn(8, 1), ...
                             'hss_ulvminnormsolve:improperRanks');

% ------------------------------------------------------- timing (information)
try
    H = mp_hss_from_generators(mp_gen_random(4096, 8192, 10, 64, 10)); rng(47); b = randn(4096, 1);
    tic; x = H \ b; tf = toc;
    tic; for r = 1:3, x = H \ b; end; ts = toc/3; %#ok<NASGU>
    fprintf('%-62s INFO  (first %.3fs, later %.4fs)\n', 'F2 4096x8192: first H\b (factor + solve), later H\b', tf, ts);
catch err
    nfail = nfail + report_error('timing run', err);
end

fprintf('\nmp_test_solver: %d failure(s)\n', nfail);
end


% ======================================================================
function f = solvecase(name, H, A, seed)
% one solve vs the dense minimum-norm solution; also checks that no tree
% level was merged.  An error counts as a failure and is reported.
try
    f = solvecase_(name, H, A, seed);
catch err
    fprintf('%-44s %-4s  (error: %s)\n', name, pf(false), err.message);
    f = 1;
end
end

function f = solvecase_(name, H, A, seed)
if isempty(A), A = full(H); end
[m, n] = size(A);
rng(seed); b = randn(m, 1);
if ~isreal(A), b = b + 1i*randn(m, 1); end
clearfactors(H);
x = H \ b;                               % factor + solve
F = H.factorcache.ulv;
xh = H \ b;                              % from the stored factors
xr = mp_minnorm_dense(A, b);
s = svd(A); kappa = s(1)/s(end);
e = mp_rel(x, xr);
bw = norm(b - A*x)/(s(1)*norm(x) + norm(b));
lv = depth_of(F);
ok = e < max(1e-11, 100*kappa*eps) && bw < 1e-14 && lv == H.levelcount && mp_rel(xh, x) == 0;
fprintf('%-44s %-11s %-9.1e %-9.1e %d/%-5d %s\n', name, sprintf('%dx%d', m, n), e, bw, ...
        lv, H.levelcount, pf(ok));
f = ~ok;
end

function f = chk(name, err, tol)
fprintf('%-62s %-4s  (%.2e)\n', name, pf(err <= tol), err);
f = ~(err <= tol);
end

function f = expect_error(name, fn, id)
try
    fn();
    fprintf('%-62s %-4s  (no error raised)\n', name, pf(false)); f = 1;
catch err
    ok = strcmp(err.identifier, id);
    fprintf('%-62s %-4s  (%s)\n', name, pf(ok), err.identifier); f = ~ok;
end
end

function f = report_error(name, err)
fprintf('%-62s %-4s  (error: %s)\n', name, pf(false), err.message);
f = 1;
end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

function d = depth_of(F)
d = 0;
while ~F.isroot
    d = d + 1;
    F = F.next;
end
end

function y = null_part(A, x)
% component of x in the null space of A (A wide, full row rank)
[Q, ~] = qr(A', 0);
y = x - Q*(Q'*x);
end

function A = nested_lowrank(m, n, minsize, k)
% copy of test_hss.nestedLowRank (hss_examples/tests/test_hss.m)
if m <= minsize || n <= minsize
    [Q1, ~] = qr(randn(m)); [Q2, ~] = qr(randn(n));
    A = Q1*[eye(m), zeros(m, n-m)]*Q2';
    return
end
m1 = ceil(m/2); m2 = m - m1; n1 = ceil(n/2); n2 = n - n1;
A11 = nested_lowrank(m1, n1, minsize, k);
A22 = nested_lowrank(m2, n2, minsize, k);
A = [A11, randn(m1, k)*randn(k, n2); randn(m2, k)*randn(k, n1), A22];
end
