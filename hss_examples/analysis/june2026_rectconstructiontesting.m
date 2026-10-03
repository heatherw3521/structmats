%% test_hss_rank_scaling.m
% Benchmarks HSS construction for a rectangular Cauchy matrix as a function
% of prescribed off-diagonal rank k (1:100), with the matrix size scaled so
% that the leaf-level blocks are always well-resolved.
%
% Three panels are produced:
%   1. Construction time  vs k
%   2. HSS storage        vs k  (dashed line = full matrix storage)
%   3. Frobenius-norm reconstruction error  vs k
%
% Requires:   hss.m / hss_constructor.m and the hss_full (or full) command.
% Usage:      test_hss_rank_scaling

clear; clc; close all;

%% ── Parameters ──────────────────────────────────────────────────────────────
k_values   = 1:100;          % ranks to test
n_rows     = 2000;           % matrix rows  (rectangular: 2000 x 1600)
n_cols     = 3000;           % matrix cols
blocksize  = 200;            % leaf block target size

fprintf('Matrix size : %d x %d\n', n_rows, n_cols);
fprintf('Blocksize   : %d\n',       blocksize);
fprintf('Ranks tested: %d – %d\n',  k_values(1), k_values(end));
fprintf('----------------------------------------------\n');

%% ── Build the Cauchy matrix once ────────────────────────────────────────────
% A(i,j) = 1 / (s(i) + t(j))   with s, t chosen so the matrix is
% well-conditioned and rectangular.
rng(42);
s = (1:n_rows)';                              % source nodes
t = (n_rows + 1 : n_rows + n_cols)';         % target nodes (disjoint → no singularity)
A_full = 1 ./ (s + t');                       % n_rows x n_cols Cauchy matrix

full_storage = numel(A_full) * 8;             % bytes (double precision)
norm_A       = norm(A_full, 'fro');

%% ── Pre-allocate result arrays ───────────────────────────────────────────────
n_k          = numel(k_values);
times        = zeros(1, n_k);
storage_hss  = zeros(1, n_k);
errors       = zeros(1, n_k);

%% ── Main loop ────────────────────────────────────────────────────────────────
fprintf('%-6s  %-12s  %-14s  %-12s\n', 'k', 'Time (s)', 'HSS storage(B)', 'Rel. error');
fprintf('%-6s  %-12s  %-14s  %-12s\n', '------','------------','------------------','----------------');
fprintf('Warming up...\n');
hss(A_full(1:400, 1:400), blocksize=blocksize, k=1);
fprintf('Done.\n\n');
for idx = 1:n_k
    k = k_values(idx);

    %── Time construction ────────────────────────────────────────────────────
    t_start = tic;
    H = hss(A_full, blocksize = blocksize, k = k);
    times(idx) = toc(t_start);

    %── Storage ──────────────────────────────────────────────────────────────
    storage_hss(idx) = hss_storage_bytes(H);

    %── Accuracy ─────────────────────────────────────────────────────────────
    A_recon        = full(H);                 % dense reconstruction
    errors(idx)    = norm(A_full - A_recon, 'fro') / norm_A;

    if mod(idx, 10) == 0 || idx == 1
        fprintf('%-6d  %-12.4f  %-14d  %-12.2e\n', ...
                k, times(idx), storage_hss(idx), errors(idx));
    end
end

%% ── Plotting ─────────────────────────────────────────────────────────────────
fig = figure('Name','HSS Rank Scaling Benchmark', ...
             'Units','normalized','Position',[0.05 0.1 0.9 0.75]);

color_hss  = [0.18 0.45 0.69];   % steel blue
color_full = [0.80 0.12 0.11];   % crimson

% ── Panel 1: Construction time ───────────────────────────────────────────────
ax1 = subplot(1,3,1);
plot(k_values, times, '-o', ...
     'Color', color_hss, 'MarkerFaceColor', color_hss, ...
     'MarkerSize', 3, 'LineWidth', 1.4);
xlabel('Off-diagonal rank  k',   'FontSize', 12);
ylabel('Construction time  (s)', 'FontSize', 12);
title('Construction Time vs Rank', 'FontSize', 13, 'FontWeight','bold');
grid on;  box on;
xlim([k_values(1), k_values(end)]);

% ── Panel 2: Storage ─────────────────────────────────────────────────────────
ax2 = subplot(1,3,2);
plot(k_values, storage_hss / 1e6, '-o', ...
     'Color', color_hss, 'MarkerFaceColor', color_hss, ...
     'MarkerSize', 3, 'LineWidth', 1.4);
hold on;
yline(full_storage / 1e6, '--', 'Full matrix', ...
      'Color', color_full, 'LineWidth', 1.8, ...
      'LabelHorizontalAlignment','left', 'FontSize', 11);
hold off;
xlabel('Off-diagonal rank  k', 'FontSize', 12);
ylabel('Storage  (MB)',         'FontSize', 12);
title('Storage vs Rank',        'FontSize', 13, 'FontWeight','bold');
legend({'HSS', 'Full matrix'}, 'Location','northwest', 'FontSize', 10);
grid on;  box on;
xlim([k_values(1), k_values(end)]);

% ── Panel 3: Reconstruction error ────────────────────────────────────────────
ax3 = subplot(1,3,3);
semilogy(k_values, errors, '-o', ...
         'Color', color_hss, 'MarkerFaceColor', color_hss, ...
         'MarkerSize', 3, 'LineWidth', 1.4);
xlabel('Off-diagonal rank  k',            'FontSize', 12);
ylabel('Relative Frobenius error',         'FontSize', 12);
title('Reconstruction Accuracy vs Rank',   'FontSize', 13, 'FontWeight','bold');
grid on;  box on;
xlim([k_values(1), k_values(end)]);

% ── Common formatting ─────────────────────────────────────────────────────────
set([ax1 ax2 ax3], 'FontSize', 11, 'TickDir', 'out');
sgtitle(sprintf('HSS Scaling  —  Cauchy matrix %d×%d,  blocksize %d', ...
                n_rows, n_cols, blocksize), ...
        'FontSize', 14, 'FontWeight', 'bold');

%saveas(fig, 'hss_rank_scaling.png');
%fprintf('\nFigure saved to hss_rank_scaling.png\n');


%% ── Helper: count HSS storage in bytes ──────────────────────────────────────
function nb = hss_storage_bytes(H)
% Recursively sum the bytes stored in all dense / low-rank fields.
    nb = 0;
    if isempty(H), return; end

    % Dense diagonal leaf block
    if ~isempty(H.D),            nb = nb + numel(H.D)            * 8; end
    % Low-rank factors and the sampled block
    if ~isempty(H.Z),            nb = nb + numel(H.Z)            * 8; end
    if ~isempty(H.Y),            nb = nb + numel(H.Y)            * 8; end
    if ~isempty(H.lrcomponent),  nb = nb + numel(H.lrcomponent)  * 8; end

    % Recurse into children
    if ~isempty(H.A11), nb = nb + hss_storage_bytes(H.A11); end
    if ~isempty(H.A22), nb = nb + hss_storage_bytes(H.A22); end
    if ~isempty(H.A12), nb = nb + hss_storage_bytes(H.A12); end
    if ~isempty(H.A21), nb = nb + hss_storage_bytes(H.A21); end
end