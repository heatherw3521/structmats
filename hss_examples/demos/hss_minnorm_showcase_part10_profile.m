%% Profiling harness for Part 10 of hss_minnorm_showcase.m
% Isolates just the large matrix-free sizes (up to 100,000 rows) under
% MATLAB's profiler, independent of the rest of the showcase file.

clear; clc; close all;
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
hssClassFile = which('hss');
if isempty(hssClassFile)
    error(['Could not locate the hss class on the path (tried adding repoRoot=%s). ' ...
        'Run this file directly so mfilename resolves correctly.'], repoRoot);
end
addpath(fullfile(repoRoot, 'hss_examples', 'lib'));   % baselines (cgne_minnorm, ...) and test helpers

%% Part 10 -- Pushing further still: genuinely large, matrix-free throughout
% Same sizes as hss_minnorm_showcase.m's Part 10; wrapped in the profiler.

ratio = 3; blocksize = 64; comptol = 1e-9;
Ms_large = [30000, 60000, 100000];
fprintf('Part 10 (profiled): matrix-free kernel, pushed larger\n');
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

function [Afun, sizeMN] = cauchyInterlacedHandle(M, ratio)
% Same as the copy in hss_minnorm_showcase.m -- interlaced Cauchy kernel
% A(i,j) = 1/(x_i - y_j) returned as a function handle so hss() never
% touches a formed dense array.
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
