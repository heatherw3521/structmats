%% june2026_conditioning_sweep.m
%
%  Compares four solvers on a rectangular HSS matrix as condition number
%  is swept from 10^1 to 10^14, with matrix size and HSS rank fixed.
%
%  Matrix construction (SVD-replacement):
%    1. Build log-kernel matrix  A_base(i,j) = log|s_i - t_j|
%       with separated grids  s in [0, 0.4],  t in [0.6, 1.0]
%       (separation guarantees low HSS off-diagonal rank)
%    2. Thin SVD:  [U, Sb, V] = svd(A_base, 'econ')
%    3. Replace singular values: sigma = logspace(0, -log10(kappa), m)
%    4. Reconstruct:  A = U * diag(sigma) * V'
%       --> exact condition number kappa, U/V inherited so HSS rank is preserved
%
%  Solvers compared:
%    (1) MATLAB backslash           A \ b
%    (2) MATLAB lsqr                lsqr(A, b, ...)
%    (3) HSS lsqr                   hss_lsqr(H, b, ...)   [HSS matvecs]
%    (4) HSS normal equations       H_B \ b  then  A'*y   [H_B = HSS(A*A')]
%
%  Four panels:
%    1. Solve time          vs log10(kappa)
%    2. lsqr iteration count vs log10(kappa)
%    3. Relative residual   vs log10(kappa)
%    4. Relative solution error vs log10(kappa)
%
%  Sanity-check figure (separate):
%    HSS Frobenius reconstruction error vs log10(kappa)
%    -- should stay flat; if it grows, rank k is insufficient

clear; clc; close all;

%% ── Parameters ───────────────────────────────────────────────────────────────
m          = 1000;          % rows  (m < n → underdetermined)
n          = 1500;         % cols
blocksize  = 200;
k_hss      = 50;           % fixed HSS off-diagonal rank

kappa_exp  = linspace(1, 14, 40);   % log10(condition number)
kappa_vals = 10 .^ kappa_exp;

lsqr_tol   = 1e-10;
lsqr_maxit = 1000;
nreps       = 3;           % timing repetitions

rng(42);

fprintf('Matrix size   : %d x %d\n', m, n);
fprintf('HSS blocksize : %d,  rank k = %d\n', blocksize, k_hss);
fprintf('kappa range   : 10^%.0f  to  10^%.0f  (%d values)\n', ...
        kappa_exp(1), kappa_exp(end), numel(kappa_vals));
fprintf('lsqr tol = %.0e,  maxit = %d,  nreps = %d\n', ...
        lsqr_tol, lsqr_maxit, nreps);
fprintf('------------------------------------------------------------------\n');

%% ── Build base log-kernel matrix and its SVD (once) ─────────────────────────
fprintf('Building log-kernel base matrix (%dx%d) and SVD...\n', m, n);

s = linspace(0.00, 0.40, m)';          % source grid  [0, 0.4]
t = linspace(0.60, 1.00, n);           % target grid  [0.6, 1.0]  (separated)
A_base = log(abs(s - t));              % m x n log-kernel matrix

[U, Sb, V] = svd(A_base, 'econ');      % U: m×m,  Sb: m×m,  V: n×m
sigma_base = diag(Sb);

fprintf('  Native condition number : %.3e\n', sigma_base(1)/sigma_base(end));
fprintf('  Numerical rank (1e-10)  : %d / %d\n\n', ...
        sum(sigma_base > 1e-10 * sigma_base(1)), m);

%% ── Warm up HSS JIT ──────────────────────────────────────────────────────────
fprintf('Warming up...\n');
sigma_w  = logspace(0, -2, m)';
A_w      = U * diag(sigma_w) * V';
H_w      = hss(A_w(1:300,1:450), blocksize=blocksize, k=k_hss);
hss_lsqr(H_w, randn(300,1), lsqr_tol, 50);
clear A_w H_w sigma_w
fprintf('Done.\n\n');

%% ── Pre-allocate ─────────────────────────────────────────────────────────────
nk = numel(kappa_vals);

time_backslash  = nan(1,nk);
time_mlsqr      = nan(1,nk);
time_hlsqr      = nan(1,nk);
time_normeq     = nan(1,nk);

iters_mlsqr     = nan(1,nk);
iters_hlsqr     = nan(1,nk);

res_backslash   = nan(1,nk);
res_mlsqr       = nan(1,nk);
res_hlsqr       = nan(1,nk);
res_normeq      = nan(1,nk);

err_backslash   = nan(1,nk);
err_mlsqr       = nan(1,nk);
err_hlsqr       = nan(1,nk);
err_normeq      = nan(1,nk);

hss_fro_err     = nan(1,nk);   % sanity check: HSS reconstruction quality
cond_check      = nan(1,nk);   % verify prescribed kappa

%% ── Main sweep ───────────────────────────────────────────────────────────────
fprintf('%-8s  %-10s  %-10s  %-10s  %-10s  %-8s  %-8s\n', ...
        'log10(k)', 't_bs', 't_mlsqr', 't_hlsqr', 't_normeq', 'it_ml', 'it_hl');
fprintf('%s\n', repmat('-',1,76));

for idx = 1:nk
    kappa = kappa_vals(idx);

    %── Construct matrix with this condition number ───────────────────────────
    sigma_new   = logspace(0, -log10(kappa), m)';
    A           = U * diag(sigma_new) * V';
    cond_check(idx) = sigma_new(1) / sigma_new(end);   % == kappa by construction

    %── Build HSS on A ────────────────────────────────────────────────────────
    H = hss(A, blocksize=blocksize, k=k_hss);

    %── HSS reconstruction sanity check ─────────────────────────────────────
    A_recon         = full(H);
    hss_fro_err(idx) = norm(A - A_recon, 'fro') / norm(A, 'fro');

    %── Build RHS: x_true in row space of A → min-norm solution is x_true ────
    x_true = V * randn(m,1);
    x_true = x_true / norm(x_true);
    b      = A * x_true;

    %── (1) MATLAB backslash ──────────────────────────────────────────────────
    t0 = tic;
    for r = 1:nreps, x1 = A \ b; end
    time_backslash(idx) = toc(t0) / nreps;
    res_backslash(idx)  = norm(A*x1 - b) / norm(b);
    err_backslash(idx)  = norm(x1 - x_true) / norm(x_true);

    %── (2) MATLAB lsqr (full matrix) ────────────────────────────────────────
    t0 = tic;
    for r = 1:nreps
        [x2, ~, ~, it2] = lsqr(A, b, lsqr_tol, lsqr_maxit);
    end
    time_mlsqr(idx)  = toc(t0) / nreps;
    iters_mlsqr(idx) = it2;
    res_mlsqr(idx)   = norm(A*x2 - b) / norm(b);
    err_mlsqr(idx)   = norm(x2 - x_true) / norm(x_true);

    %── (3) HSS lsqr ─────────────────────────────────────────────────────────
    t0 = tic;
    for r = 1:nreps
        [x3, ~, ~, it3] = hss_lsqr(H, b, lsqr_tol, lsqr_maxit);
    end
    time_hlsqr(idx)  = toc(t0) / nreps;
    iters_hlsqr(idx) = it3;
    res_hlsqr(idx)   = norm(A*x3 - b) / norm(b);
    err_hlsqr(idx)   = norm(x3 - x_true) / norm(x_true);

    %── (4) HSS normal equations ──────────────────────────────────────────────
    B   = A * A';
    HB  = hss(B, blocksize=blocksize);
    t0  = tic;
    for r = 1:nreps
        y  = HB \ b;
        x4 = A' * y;
    end
    time_normeq(idx) = toc(t0) / nreps;
    res_normeq(idx)  = norm(A*x4 - b) / norm(b);
    err_normeq(idx)  = norm(x4 - x_true) / norm(x_true);

    if mod(idx,5)==0 || idx==1
        fprintf('%-8.2f  %-10.3f  %-10.3f  %-10.3f  %-10.3f  %-8d  %-8d\n', ...
                kappa_exp(idx), ...
                time_backslash(idx), time_mlsqr(idx), ...
                time_hlsqr(idx),     time_normeq(idx), ...
                iters_mlsqr(idx),    iters_hlsqr(idx));
    end

    clear A B HB x1 x2 x3 x4 y H
end

%% ── Save results ─────────────────────────────────────────────────────────────
save('june2026_conditioning_sweep.mat', ...
     'kappa_exp','kappa_vals','m','n','blocksize','k_hss', ...
     'time_backslash','time_mlsqr','time_hlsqr','time_normeq', ...
     'iters_mlsqr','iters_hlsqr', ...
     'res_backslash','res_mlsqr','res_hlsqr','res_normeq', ...
     'err_backslash','err_mlsqr','err_hlsqr','err_normeq', ...
     'hss_fro_err','cond_check');
fprintf('\nResults saved.\n\n');

%% ── Plotting ─────────────────────────────────────────────────────────────────

% ── Appearance (matches may2026 style) ───────────────────────────────────────
FONT      = 'Times New Roman';
FONT_SZ   = 11;
LABEL_SZ  = 12;
TITLE_SZ  = 12;
LEGEND_SZ = 10;
LW        = 1.6;
MS        = 5;

% Wong (2011) colour-blind-safe palette
C_bs  = [0   114 178]/255;   % blue       – backslash
C_ml  = [230 159   0]/255;   % orange     – MATLAB lsqr
C_hl  = [0   158 115]/255;   % green      – HSS lsqr
C_ne  = [213  94   0]/255;   % vermilion  – HSS normeq

colors  = {C_bs, C_ml, C_hl, C_ne};
styles  = {'-o', '--s', '-.^', ':d'};
labels  = {'MATLAB $\backslash$', 'MATLAB lsqr', 'HSS lsqr', 'HSS norm-eq'};

%% ── Figure 1: Main benchmark (2×2) ──────────────────────────────────────────
fig1 = figure('Name','Conditioning Sweep', ...
              'Units','centimeters','Position',[2 2 18 14]);

time_data = {time_backslash, time_mlsqr, time_hlsqr, time_normeq};
iter_data = {iters_mlsqr, iters_hlsqr};
res_data  = {res_backslash,  res_mlsqr,  res_hlsqr,  res_normeq};
err_data  = {err_backslash,  err_mlsqr,  err_hlsqr,  err_normeq};

% Panel 1: Solve time
ax1 = subplot(2,2,1);
for s = 1:4
    semilogy(kappa_exp, time_data{s}, styles{s}, ...
             'Color',colors{s}, 'LineWidth',LW, ...
             'MarkerSize',MS, 'MarkerFaceColor',colors{s}, ...
             'DisplayName',labels{s});
    hold on;
end
hold off;
ylabel('Avg solve time (s)', 'FontSize',LABEL_SZ, 'FontName',FONT);
title('Solve Time', 'FontSize',TITLE_SZ, 'FontName',FONT, 'FontWeight','bold');
legend('Location','northwest','FontSize',LEGEND_SZ,'FontName',FONT, ...
       'Interpreter','latex','Box','on');
style_ax(ax1, FONT_SZ, FONT);

% Panel 2: Iteration count (lsqr methods only)
ax2 = subplot(2,2,2);
iter_styles = {'--s', '-.^'};
iter_colors = {C_ml, C_hl};
iter_labels = {'MATLAB lsqr', 'HSS lsqr'};
for s = 1:2
    plot(kappa_exp, iter_data{s}, iter_styles{s}, ...
         'Color',iter_colors{s}, 'LineWidth',LW, ...
         'MarkerSize',MS, 'MarkerFaceColor',iter_colors{s}, ...
         'DisplayName',iter_labels{s});
    hold on;
end
yline(lsqr_maxit, ':k', 'maxit', 'LineWidth',1.1, ...
      'LabelHorizontalAlignment','left', 'FontSize',FONT_SZ-1);
hold off;
ylabel('Iterations to converge', 'FontSize',LABEL_SZ, 'FontName',FONT);
title('lsqr Iteration Count', 'FontSize',TITLE_SZ, 'FontName',FONT, 'FontWeight','bold');
legend('Location','northwest','FontSize',LEGEND_SZ,'FontName',FONT,'Box','on');
style_ax(ax2, FONT_SZ, FONT);

% Panel 3: Relative residual
ax3 = subplot(2,2,3);
for s = 1:4
    semilogy(kappa_exp, res_data{s}, styles{s}, ...
             'Color',colors{s}, 'LineWidth',LW, ...
             'MarkerSize',MS, 'MarkerFaceColor',colors{s}, ...
             'DisplayName',labels{s});
    hold on;
end
semilogy(kappa_exp, eps * kappa_vals, '--', ...
         'Color',[0.6 0.6 0.6], 'LineWidth',1.1, ...
         'DisplayName','$\epsilon\kappa$');
hold off;
xlabel('$\log_{10}(\kappa)$', 'Interpreter','latex', ...
       'FontSize',LABEL_SZ, 'FontName',FONT);
ylabel('$\|Ax-b\|/\|b\|$', 'Interpreter','latex', ...
       'FontSize',LABEL_SZ, 'FontName',FONT);
title('Relative Residual', 'FontSize',TITLE_SZ, 'FontName',FONT, 'FontWeight','bold');
legend('Location','northwest','FontSize',LEGEND_SZ,'FontName',FONT, ...
       'Interpreter','latex','Box','on');
style_ax(ax3, FONT_SZ, FONT);

% Panel 4: Relative solution error
ax4 = subplot(2,2,4);
for s = 1:4
    semilogy(kappa_exp, err_data{s}, styles{s}, ...
             'Color',colors{s}, 'LineWidth',LW, ...
             'MarkerSize',MS, 'MarkerFaceColor',colors{s}, ...
             'DisplayName',labels{s});
    hold on;
end
semilogy(kappa_exp, eps * kappa_vals,   '--', 'Color',[0.6 0.6 0.6], ...
         'LineWidth',1.1, 'DisplayName','$\epsilon\kappa$');
semilogy(kappa_exp, eps * kappa_vals.^2, ':', 'Color',[0.6 0.6 0.6], ...
         'LineWidth',1.1, 'DisplayName','$\epsilon\kappa^2$');
hold off;
xlabel('$\log_{10}(\kappa)$', 'Interpreter','latex', ...
       'FontSize',LABEL_SZ, 'FontName',FONT);
ylabel('$\|x-x_{\mathrm{true}}\|/\|x_{\mathrm{true}}\|$', ...
       'Interpreter','latex', 'FontSize',LABEL_SZ, 'FontName',FONT);
title('Relative Solution Error', 'FontSize',TITLE_SZ, 'FontName',FONT, 'FontWeight','bold');
legend('Location','northwest','FontSize',LEGEND_SZ,'FontName',FONT, ...
       'Interpreter','latex','Box','on');
style_ax(ax4, FONT_SZ, FONT);

sgtitle(sprintf('Solver comparison — log-kernel %d\\times%d,  blocksize %d,  HSS rank k=%d', ...
                m, n, blocksize, k_hss), ...
        'FontSize',TITLE_SZ+1, 'FontName',FONT, 'FontWeight','bold');

for ax = [ax1 ax2 ax3 ax4]
    xlim(ax, [kappa_exp(1) kappa_exp(end)]);
end

%export_fig_pub(fig1, 'june2026_conditioning_sweep_main');

%% ── Figure 2: HSS sanity check ───────────────────────────────────────────────
fig2 = figure('Name','HSS Reconstruction Sanity Check', ...
              'Units','centimeters','Position',[22 2 10 7]);

semilogy(kappa_exp, hss_fro_err, '-o', ...
         'Color',[0 158 115]/255, 'LineWidth',LW, ...
         'MarkerSize',MS, 'MarkerFaceColor',[0 158 115]/255);
xlabel('$\log_{10}(\kappa)$', 'Interpreter','latex', ...
       'FontSize',LABEL_SZ, 'FontName',FONT);
ylabel('$\|A - A_{\mathrm{HSS}}\|_F / \|A\|_F$', ...
       'Interpreter','latex', 'FontSize',LABEL_SZ, 'FontName',FONT);
title({'HSS Reconstruction Error vs Condition Number'; ...
       sprintf('(if flat: rank k=%d is adequate for all \\kappa)', k_hss)}, ...
      'FontSize',TITLE_SZ, 'FontName',FONT, 'FontWeight','bold');
xlim([kappa_exp(1) kappa_exp(end)]);
style_ax(gca, FONT_SZ, FONT);

%export_fig_pub(fig2, 'june2026_conditioning_sweep_sanity');

%% ── Console summary ──────────────────────────────────────────────────────────
fprintf('\n── Summary ──────────────────────────────────────────────────────\n');

thresh_err = 1e-2;
% for [data, name] = [err_backslash; {"backslash"}; ...
%                     err_mlsqr;    {"MATLAB lsqr"}; ...
%                     err_hlsqr;    {"HSS lsqr"}; ...
%                     err_normeq;   {"HSS normeq"}]
%     % (plain loop is cleaner in MATLAB)
% end
solver_names = {'backslash', 'MATLAB lsqr', 'HSS lsqr', 'HSS normeq'};
solver_errs  = {err_backslash, err_mlsqr, err_hlsqr, err_normeq};
for s = 1:4
    bad = solver_errs{s} > thresh_err;
    if any(bad)
        fprintf('%s: solution error > %.0e at log10(kappa) >= %.1f\n', ...
                solver_names{s}, thresh_err, kappa_exp(find(bad,1)));
    else
        fprintf('%s: solution error stayed below %.0e across all kappa.\n', ...
                solver_names{s}, thresh_err);
    end
end

if all(hss_fro_err < 1e-4)
    fprintf('\nHSS sanity: reconstruction error < 1e-4 for all kappa. Rank k=%d is adequate.\n', k_hss);
else
    first_bad = kappa_exp(find(hss_fro_err >= 1e-4, 1));
    fprintf('\nHSS sanity WARNING: reconstruction error >= 1e-4 from log10(kappa) = %.1f.\n', first_bad);
    fprintf('  Consider increasing k_hss above %d.\n', k_hss);
end
fprintf('─────────────────────────────────────────────────────────────────\n');

%% ── Local helpers ────────────────────────────────────────────────────────────

function [x, flag, relres, iter] = hss_lsqr(H, b, tol, maxit)
    fun        = @(x, mode) lsqr_wrapper(x, mode, H);
    [x, flag, relres, iter] = lsqr(fun, b, tol, maxit);
end

function y = lsqr_wrapper(x, mode, H)
    if strcmp(mode, 'notransp')
        y = H * x;
    else
        y = H.' * x;
    end
end

function style_ax(ax, fs, font)
    set(ax, 'FontSize',fs, 'FontName',font, ...
            'XColor','black', 'YColor','black', ...
            'Box','on', 'TickDir','out', ...
            'XMinorTick','on', 'YMinorTick','on', ...
            'LineWidth',0.8, 'GridAlpha',0.25, ...
            'GridLineStyle','--');
    grid(ax, 'on');
end

function export_fig_pub(fig, fname)
    exportgraphics(fig, [fname '.pdf'], ...
        'ContentType','vector','BackgroundColor','white');
    exportgraphics(fig, [fname '.png'], ...
        'Resolution',300,'BackgroundColor','white');
    fprintf('  Saved  %s.{pdf,png}\n', fname);
end