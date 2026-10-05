%% HSS Underdetermined (Minimum-Norm) Solver Showcase
%
% Demonstrates the rectangular/minimum-norm ULV solve path -- H \ b for a
% wide (m < n) hss matrix H, which dispatches to
% @hss/private/hss_ulvminnormsolve.m -- recovering
%
%       x* = argmin ||x||_2  s.t.  Hx = b
%
% starting from a tiny 2-3 level example and building up to large,
% depth-controlled matrices. At every scale we check TWO things against
% the true minimum-norm solution (via MATLAB's lsqminnorm on the dense
% matrix):
%   1. constraint satisfaction:      ||Hx - b|| / ||b||
%   2. minimum-norm-ness:            ||x - x_minnorm|| / ||x_minnorm||
% and we time the HSS solve throughout.
%
% Test matrices are built by rectsampler() (local function below): start
% from an exactly rank-k matrix (so the HSS off-diagonal compression is
% exact), then perturb every LEAF's dense diagonal block
% with independent noise via diagmodify() so the local/diagonal blocks
% are well-conditioned (full rank) -- the one precondition this ULV-style
% solver needs, same as any HSS direct solver. Parts 1-5 use that
% construction, with an exact k always handed to both rectsampler() and
% (implicitly) to hss(). Part 6 onward switches to a genuine well-separated
% kernel matrix instead -- not built FOR the HSS format, off-diagonal rank
% discovered by tolerance rather than known in advance -- to test the
% min-norm solver against something closer to a real use case.

clear; clc; close all;
% @hss and +hssutil live at the repo root, two levels up from this file.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
hssClassFile = which('hss');
if isempty(hssClassFile)
    error(['Could not locate the hss class on the path (tried adding repoRoot=%s). ' ...
        'Run this file directly (not pasted/section-by-section) so mfilename resolves correctly.'], repoRoot);
end
addpath(fullfile(repoRoot, 'hss_examples', 'lib'));   % baselines (cgne_minnorm, ...) and test helpers

set(groot,'defaultAxesFontSize',12);
set(groot,'defaultAxesLineWidth',1.1);
set(groot,'defaultLineLineWidth',1.8);
set(groot,'defaultFigureColor','w');

TOL = 1e-8;   % pass/fail threshold for both checks below

%% Part 1 -- Simple case: 40 x 60, depth 2-3
% Smallest case worth checking by hand.

rng(1);
M = 40; N = 60; blocksize = 15; k = 6;

[A, H] = rectsampler(M, N, k, blocksize);
fprintf('Part 1: %d x %d, blocksize=%d, k=%d  ->  levelcount = %d\n', ...
    M, N, blocksize, k, H.levelcount);

b = randn(M, 1);

tic;
x = H \ b;
t_hss = toc;

tic;
xref = lsqminnorm(A, b);
t_dense = toc;

report_result('Part 1 (40x60, depth ~2-3)', A, b, x, xref, t_hss, t_dense, TOL);

%% Part 2 -- Building up: a few more depths, still small enough to be instant
% Same construction, deeper trees each time.

rng(2);
configs2 = { ...
    struct('M',  80, 'N', 120, 'blocksize', 15, 'k', 3), ...  % depth ~3
    struct('M', 160, 'N', 240, 'blocksize', 15, 'k', 3), ...  % depth ~4
    struct('M', 320, 'N', 480, 'blocksize', 15, 'k', 3)  ...  % depth ~5
};

for i = 1:numel(configs2)
    c = configs2{i};
    [A, H] = rectsampler(c.M, c.N, c.k, c.blocksize);
    b = randn(c.M, 1);

    tic; x = H \ b; t_hss = toc;
    tic; xref = lsqminnorm(A, b); t_dense = toc;

    label = sprintf('Part 2.%d (%dx%d, depth %d)', i, c.M, c.N, H.levelcount);
    report_result(label, A, b, x, xref, t_hss, t_dense, TOL);
end

%% Part 3 -- Large-scale, depth-controlled sweep
% Fix blocksize/k, grow N: depth-controlled correctness + timing sweep.

blocksize = 16;
k = 3;
Ns = [128, 256, 512, 1024, 2048, 4096];   % dense lsqminnorm stays cheap here

nN = numel(Ns);
levelcounts   = zeros(nN,1);
t_hss_solve   = zeros(nN,1);
t_dense_solve = zeros(nN,1);
residuals     = zeros(nN,1);
minnormerrs   = zeros(nN,1);
passed        = false(nN,1);

for i = 1:nN
    N = Ns(i); M = round(N/2);
    rng(10 + i);

    [A, H] = rectsampler(M, N, k, blocksize);
    b = randn(M, 1);

    tic; x = H \ b; t_hss_solve(i) = toc;
    tic; xref = lsqminnorm(A, b); t_dense_solve(i) = toc;

    levelcounts(i)  = H.levelcount;
    residuals(i)    = norm(A*x - b) / norm(b);
    minnormerrs(i)  = norm(x - xref) / norm(xref);
    passed(i)       = residuals(i) < TOL && minnormerrs(i) < TOL;
end

T = table(Ns(:), levelcounts, t_hss_solve, t_dense_solve, residuals, minnormerrs, passed, ...
    'VariableNames', {'N','levelcount','HSS_solve_s','Dense_lsqminnorm_s','Residual','MinNormError','Pass'});
disp(T);

fprintf('\nPart 3: %d / %d depth-controlled cases matched the true minimum-norm solution.\n', ...
    sum(passed), nN);

%% Part 4 -- Pushing further: depths where dense lsqminnorm starts to hurt
% Same construction, larger still: HSS stays fast, dense gets slow.

Ns_big = [8192, 16384];
nBig = numel(Ns_big);
levelcounts_big = zeros(nBig,1);
t_hss_big       = zeros(nBig,1);
t_dense_big     = zeros(nBig,1);
residuals_big   = zeros(nBig,1);
minnormerrs_big = zeros(nBig,1);
passed_big      = false(nBig,1);

for i = 1:nBig
    N = Ns_big(i); M = round(N/2);
    rng(50 + i);

    [A, H] = rectsampler(M, N, k, blocksize);
    b = randn(M, 1);

    tic; x = H \ b; t_hss_big(i) = toc;
    tic; xref = lsqminnorm(A, b); t_dense_big(i) = toc;

    levelcounts_big(i) = H.levelcount;
    residuals_big(i)   = norm(A*x - b) / norm(b);
    minnormerrs_big(i) = norm(x - xref) / norm(xref);
    passed_big(i)      = residuals_big(i) < TOL && minnormerrs_big(i) < TOL;

    fprintf('N=%6d  levelcount=%2d  HSS solve=%.3fs  dense lsqminnorm=%.3fs  residual=%.2e  errvsref=%.2e  %s\n', ...
        N, levelcounts_big(i), t_hss_big(i), t_dense_big(i), residuals_big(i), minnormerrs_big(i), ...
        pass_str(passed_big(i)));
end

%% Part 5 -- Timing plot
% Solve time vs N (Parts 3-4), HSS vs dense, log-log.

allN     = [Ns(:); Ns_big(:)];
allTHss  = [t_hss_solve; t_hss_big];
allTDense = [t_dense_solve; t_dense_big];
[allN, order] = sort(allN);
allTHss = allTHss(order);
allTDense = allTDense(order);

figure('Position', [100 100 640 480]);
loglog(allN, allTHss, '-o', 'DisplayName', 'HSS min-norm solve');
hold on;
loglog(allN, allTDense, '-s', 'DisplayName', 'Dense lsqminnorm');
grid on; box on;
xlabel('N (columns)');
ylabel('Solve time (s)');
title('HSS Minimum-Norm Solve: Timing vs. Problem Size');
legend('Location', 'northwest');

fprintf('\nAll parts complete. %d / %d total depth-controlled cases (Parts 3-4) matched the true minimum-norm solution.\n', ...
    sum(passed) + sum(passed_big), nN + nBig);

%% Part 6 -- A natural matrix, not a synthetic one
% Interlaced Cauchy kernel A(i,j)=1/(x_i-y_j), well-conditioned (cond~2),
% off-diagonal rank discovered by tolerance rather than handed to hss().

blocksize = 32; comptol = 1e-9;
Ms = [128, 512, 2048];
ratio = 3;   % y-points interlaced between each adjacent pair of x-points

fprintf('Part 6: interlaced Cauchy kernel, blocksize=%d, tol=%.0e -- natural, not synthetic\n', ...
    blocksize, comptol);
fprintf('%-8s %-8s %-10s %-6s %-10s %-10s %-10s %-10s %-10s %-10s\n', ...
    'M','N','cond(A)','levels','comprerr','residual','errvsref','build(s)','HSS(s)','dense(s)');

for i = 1:numel(Ms)
    M = Ms(i);
    A = cauchyInterlaced(M, ratio);
    N = size(A, 2);
    condA = cond(A);

    tic; H = hss(A, blocksize = blocksize, tol = comptol); tbuild = toc;
    comprerr = norm(full(H) - A, 'fro') / norm(A, 'fro');

    rng(i);
    b = randn(M, 1);
    tic; x = H \ b; t_hss = toc;
    tic; xref = lsqminnorm(A, b); t_dense = toc;
    res = norm(A*x - b) / norm(b);
    err = norm(x - xref) / norm(xref);

    fprintf('%-8d %-8d %-10.2e %-6d %-10.2e %-10.2e %-10.2e %-10.2f %-10.4f %-10.4f  %s\n', ...
        M, N, condA, H.levelcount, comprerr, res, err, tbuild, t_hss, t_dense, ...
        pass_str(res < TOL && err < TOL));
end

%% Part 7 -- Pushing the same natural matrix to real depth
% b = H*(H'*z) makes x's exact answer known, so the solve check needs
% no dense array; compression is checked separately against the true A.

Ms_deep = [2048, 4096, 8192];

fprintf('\nPart 7: interlaced Cauchy kernel\n');
fprintf('%-8s %-8s %-6s %-10s %-10s %-10s %-10s %-10s\n', ...
    'M','N','levels','build(s)','HSS(s)','comprerr','hss_res','errvstrue');

for i = 1:numel(Ms_deep)
    M = Ms_deep(i);
    A = cauchyInterlaced(M, ratio);
    N = size(A, 2);

    tic; H = hss(A, blocksize = blocksize, tol = comptol); tbuild = toc;
    comprerr = norm(full(H) - A, 'fro') / norm(A, 'fro');   % vs the TRUE A

    rng(100 + i);
    z = randn(M, 1);
    xtrue = H' * z;      % in the row space of H by construction
    b = H * xtrue;
    tic; x = H \ b; t_hss = toc;
    res = norm(H*x - b) / norm(b);           % against H itself
    errvstrue = norm(x - xtrue) / norm(xtrue); % against the exact min-norm solution of H's own system

    fprintf('%-8d %-8d %-6d %-10.2f %-10.4f %-10.2e %-10.2e %-10.2e\n', ...
        M, N, H.levelcount, tbuild, t_hss, comprerr, res, errvstrue);
end

%% Part 8 -- A genuinely different solver, for comparison: matvec-only CG
% cgne_minnorm (CG on H*H', matvec-only, never full(H)) vs. direct ULV.
% Fair comparison because this kernel is deliberately well-conditioned.

fprintf('\nPart 8: same kernel -- direct ULV solve (H\\b) vs. matvec-only CG (cgne_minnorm), accuracy AND timing\n');
fprintf('%-8s %-8s %-6s %-10s %-10s %-10s %-10s\n', ...
    'M','N','levels','ULV(s)','pcg(s)','pcg/ULV','agree');

Ms_pcg = [512, 2048, 4096, 8192];

for i = 1:numel(Ms_pcg)
    M = Ms_pcg(i);
    A = cauchyInterlaced(M, ratio);
    N = size(A, 2);
    H = hss(A, blocksize = blocksize, tol = comptol);

    rng(200 + i);
    b = randn(M, 1);

    tic; x_ulv = H \ b; t_ulv = toc;
    tic; x_pcg = cgne_minnorm(H, b); t_pcg = toc;
    agree = norm(x_ulv - x_pcg) / norm(x_ulv);

    fprintf('%-8d %-8d %-6d %-10.4f %-10.4f %-10.2f %-10.2e\n', ...
        M, N, H.levelcount, t_ulv, t_pcg, t_pcg / t_ulv, agree);
end

%% Part 9 -- The dense array never exists at all, not even during construction
% hss(Afun, sizeA=[M,N], ...) evaluates only the entries it samples.

M = 20000; ratio = 3; blocksize = 64; comptol = 1e-9;
[Afun, sizeMN] = cauchyInterlacedHandle(M, ratio);

fprintf('\nPart 9: same kernel, but the M x N array is never formed\n');
fprintf('sizeA=%s (a dense array this size would be ~%.1fGB; never allocated)\n', ...
    mat2str(sizeMN), prod(sizeMN)*8/1e9);

tic; H = hss(Afun, blocksize = blocksize, tol = comptol, sizeA = sizeMN); tbuild = toc;

rng(1);
z = randn(M, 1);
xtrue = H' * z;
b = H * xtrue;
tic; x_ulv = H \ b; t_ulv = toc;
tic; x_pcg = cgne_minnorm(H, b); t_pcg = toc;

fprintf('%-8s %-6s %-10s %-10s %-10s %-10s %-10s\n', ...
    'levels','build(s)','ULV(s)','pcg(s)','ulv_res','ulv_err','agree');
fprintf('%-8d %-10.2f %-10.4f %-10.4f %-10.2e %-10.2e %-10.2e\n', ...
    H.levelcount, tbuild, t_ulv, t_pcg, norm(H*x_ulv-b)/norm(b), ...
    norm(x_ulv-xtrue)/norm(xtrue), norm(x_ulv-x_pcg)/norm(x_ulv));

%% Part 10 -- Pushing further still: genuinely large, matrix-free throughout
% Same as Part 9, several sizes larger yet (up to 100,000 rows).

Ms_large = [30000, 60000, 100000];

fprintf('\nPart 10: matrix-free kernel, pushed larger\n');
fprintf('%-8s %-6s %-10s %-10s %-10s %-10s %-10s\n', ...
    'M','levels','build(s)','ULV(s)','pcg(s)','ulv_res','ulv_err');

for i = 1:numel(Ms_large)
    M = Ms_large(i);
    [Afun, sizeMN] = cauchyInterlacedHandle(M, ratio);

    tic; H = hss(Afun, blocksize = blocksize, tol = comptol, sizeA = sizeMN); tbuild = toc;

    rng(300 + i);
    z = randn(M, 1);
    xtrue = H' * z;
    b = H * xtrue;
    tic; x_ulv = H \ b; t_ulv = toc;
    tic; x_pcg = cgne_minnorm(H, b); t_pcg = toc;

    fprintf('%-8d %-6d %-10.2f %-10.4f %-10.4f %-10.2e %-10.2e\n', ...
        M, H.levelcount, tbuild, t_ulv, t_pcg, norm(H*x_ulv-b)/norm(b), norm(x_ulv-xtrue)/norm(xtrue));
end


%% ======================= Local functions =======================

function report_result(label, A, b, x, xref, t_hss, t_dense, tol)
res = norm(A*x - b) / norm(b);
err = norm(x - xref) / norm(xref);
ok = res < tol && err < tol;
fprintf(['%-28s  residual=%.2e  errvsref=%.2e  ||x||=%.4f  ||x_minnorm||=%.4f  ' ...
    'HSS=%.4fs  dense=%.4fs  %s\n'], ...
    label, res, err, norm(x), norm(xref), t_hss, t_dense, pass_str(ok));
end

function s = pass_str(ok)
if ok
    s = 'PASS';
else
    s = 'FAIL';
end
end

function [Z, HZ] = rectsampler(M, N, k, blocksize)
% Build an exactly rank-k wide matrix (so the HSS off-diagonal
% compression is exact -- zero
% reconstruction error, independent of tolerance/rank-detection), then
% perturb each leaf's OWN dense diagonal block with independent noise so
% the local blocks are well-conditioned (full rank). The off-diagonal
% low-rank factors are untouched by diagmodify, so full(HZ) below is
% exact relative to what HZ itself represents.
Z = rand(M, k) * rand(k, N);
HZ = hss(Z, blocksize = blocksize);
HZ = diagmodify(HZ);
Z = full(HZ);
end

function H = diagmodify(H)
if H.isleaf
    H.D = H.D + rand(size(H.D));
else
    H.A11 = diagmodify(H.A11);
    H.A22 = diagmodify(H.A22);
end
end

function A = cauchyInterlaced(M, ratio)
% Rectangular Cauchy kernel A(i,j) = 1/(x_i - y_j), with x and y
% INTERLACED along the same line rather than placed on two disjoint
% intervals. M x-points at integer positions 1..M; `ratio` y-points
% inserted strictly between every adjacent pair of x-points (so
% N = (M-1)*ratio), each offset a fraction of the gap to its neighbors so
% no y ever coincides exactly with an x (which would give a literal
% division by zero).
%
% Interlacing gives every point near neighbors (handled densely by the
% leaves) as well as far ones (compressed off-diagonal blocks), so cond(A)
% stays close to 1 at every size; points on two disjoint intervals would
% make A severely ill-conditioned (cond ~ 1e18).
x = (1:M)';
y = zeros((M-1)*ratio, 1);
idx = 1;
for i = 1:M-1
    for r = 1:ratio
        y(idx) = x(i) + r/(ratio+1) * (x(i+1) - x(i));
        idx = idx + 1;
    end
end
A = 1 ./ (x - y.');
end

function [Afun, sizeMN] = cauchyInterlacedHandle(M, ratio)
% Same interlaced Cauchy kernel as cauchyInterlaced(), but returns a
% function handle Afun(I,J) that evaluates entries on demand, for use as
% hss(Afun, blocksize=..., sizeA=sizeMN, tol=...). sizeA is required with a
% handle. The dense array is never formed; only the blocks the constructor
% samples are evaluated. See Part 9.
N = (M-1)*ratio;
x = (1:M)';
y = zeros(N, 1);
idx = 1;
for i = 1:M-1
    for r = 1:ratio
        y(idx) = x(i) + r/(ratio+1) * (x(i+1) - x(i));
        idx = idx + 1;
    end
end
Afun = @(I,J) 1 ./ (x(I) - y(J).');
sizeMN = [M, N];
end
