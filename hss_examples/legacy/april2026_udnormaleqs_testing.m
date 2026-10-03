%% Rectangular HSS underdetermined solve test via normal equations (A*A^T)
clear; clc;

% -----------------------------------------------------------------------
% Parameters
% -----------------------------------------------------------------------
Ms            = round(logspace(3, 4, 8));
aspect_ratios = [1.5, 2.0, 2.5, 3.0];
k             = 45;
nreps         = 3;
max_full_size = 12000;

results = struct();

for ar_idx = 1:length(aspect_ratios)
    aspect_ratio = aspect_ratios(ar_idx);
    N_sizes      = round(Ms * aspect_ratio);

    fprintf('\n########################################\n');
    fprintf('ASPECT RATIO: %.1f:1  (N > M, underdetermined)\n', aspect_ratio);
    fprintf('########################################\n');

    % pre-allocate
    hss_build_time      = nan(size(Ms));
    matlab_lsqr_time    = nan(size(Ms));
    hss_normeq_time     = nan(size(Ms));

    matlab_residual     = nan(size(Ms));   % ||A x_ml - b||
    hss_normeq_residual = nan(size(Ms));   % ||A x_hss - b||
    normeq_error        = nan(size(Ms));   % ||x_hss - x_ml||  (solution error)

    matlab_xnorm        = nan(size(Ms));   % ||x_ml||
    hss_xnorm           = nan(size(Ms));   % ||x_hss||

    full_storage        = nan(size(Ms));
    gram_storage        = nan(size(Ms));
    hss_storage         = nan(size(Ms));

    for idx = 1:length(Ms)
        M = Ms(idx);
        N = N_sizes(idx);

        fprintf('\n====================================\n');
        fprintf('Size: M = %d, N = %d  (%.1fx underdetermined)\n', M, N, N/M);
        fprintf('====================================\n');

        can_form_full = (M * N < max_full_size^2);

        if ~can_form_full
            fprintf('Skipping: A (%d x %d) exceeds full-matrix size limit.\n', M, N);
            continue
        end

        % ---- Build A ----------------------------------------------------
        fprintf('Constructing full A (%d x %d)...\n', M, N);
        s = (1:M)' / M;
        t = (1:N)  / N;
        A = 1 ./ (s + t);

        full_info         = whos('A');
        full_storage(idx) = full_info.bytes / 1e6;

        % ---- Form B = A*A' ----------------------------------------------
        fprintf('Forming B = A * A''  (%d x %d)...\n', M, M);
        B = A * A';

        gram_info         = whos('B');
        gram_storage(idx) = gram_info.bytes / 1e6;

        % ---- Build HSS(B) -----------------------------------------------
        blocksize = min(max(200, round(0.5*sqrt(M))), 2000);
        fprintf('Building HSS(B) with blocksize %d...\n', blocksize);

        t0 = tic;
        HB = hss(B, blocksize = blocksize);
        hss_build_time(idx) = toc(t0);

        hss_info          = whos('HB');
        hss_storage(idx)  = hss_info.bytes / 1e6;
        fprintf('HSS build time: %.3f s\n', hss_build_time(idx));

        % ---- Test problem ------------------------------------------------
        x_true = randn(N, 1);
        x_true = x_true / norm(x_true);
        b      = A * x_true;

        % ================================================================
        % Method 1: MATLAB lsqr on full A  (baseline)
        % ================================================================
        fprintf('Testing MATLAB lsqr...\n');
        tic;
        for rep = 1:nreps
            x_ml = lsqr(A, b, 1e-10, 2000);
        end
        matlab_lsqr_time(idx) = toc / nreps;
        matlab_residual(idx)  = norm(A * x_ml - b);
        matlab_xnorm(idx)     = norm(x_ml);

        fprintf('Avg MATLAB lsqr:   %.4f s  |  ||Ax-b|| = %.2e  |  ||x|| = %.4e\n', ...
                matlab_lsqr_time(idx), matlab_residual(idx), matlab_xnorm(idx));

        % ================================================================
        % Method 2: HSS normal equations
        %   Solve  (A*A') y = b  via HSS backslash, then x = A'*y
        % ================================================================
        fprintf('Testing HSS normal-eq solve...\n');
        tic;
        for rep = 1:nreps
            y_hss = HB \ b;
            x_hss = A' * y_hss;
        end
        hss_normeq_time(idx)     = toc / nreps;
        hss_normeq_residual(idx) = norm(A * x_hss - b);
        hss_xnorm(idx)           = norm(x_hss);
        normeq_error(idx)        = norm(x_hss - x_ml);   % vs MATLAB baseline

        fprintf('Avg HSS normal-eq: %.4f s  |  ||Ax-b|| = %.2e  |  ||x|| = %.4e\n', ...
                hss_normeq_time(idx), hss_normeq_residual(idx), hss_xnorm(idx));
        fprintf('  Solution error vs MATLAB:  ||x_hss - x_ml|| = %.2e\n', normeq_error(idx));
        fprintf('  Min-norm check:            ||x_ml|| = %.4e   ||x_hss|| = %.4e\n', ...
                matlab_xnorm(idx), hss_xnorm(idx));

        clear B HB A y_hss x_hss x_ml
    end

    % -------------------------------------------------------------------
    % Store results
    % -------------------------------------------------------------------
    results(ar_idx).aspect_ratio        = aspect_ratio;
    results(ar_idx).M                   = Ms(:);
    results(ar_idx).N                   = N_sizes(:);
    results(ar_idx).hss_build_time      = hss_build_time(:);
    results(ar_idx).matlab_lsqr_time    = matlab_lsqr_time(:);
    results(ar_idx).hss_normeq_time     = hss_normeq_time(:);
    results(ar_idx).matlab_residual     = matlab_residual(:);
    results(ar_idx).hss_normeq_residual = hss_normeq_residual(:);
    results(ar_idx).normeq_error        = normeq_error(:);
    results(ar_idx).matlab_xnorm        = matlab_xnorm(:);
    results(ar_idx).hss_xnorm           = hss_xnorm(:);
    results(ar_idx).full_storage        = full_storage(:);
    results(ar_idx).gram_storage        = gram_storage(:);
    results(ar_idx).hss_storage         = hss_storage(:);

    results(ar_idx).table = table( ...
        Ms(:), N_sizes(:), ...
        hss_build_time(:), ...
        matlab_lsqr_time(:),  hss_normeq_time(:), ...
        matlab_residual(:),   hss_normeq_residual(:), normeq_error(:), ...
        matlab_xnorm(:),      hss_xnorm(:), ...
        'VariableNames', { ...
            'M', 'N', ...
            'HSS_build_s', ...
            'MATLAB_lsqr_s',  'HSS_normeq_s', ...
            'MATLAB_residual', 'HSS_normeq_residual', 'HSS_solution_error', ...
            'MATLAB_xnorm',    'HSS_xnorm'});

    fname = sprintf('hss_underdetermined_AR%.1f.txt', aspect_ratio);
    writetable(results(ar_idx).table, fname, 'Delimiter', '\t');
    fprintf('\nSaved results to %s\n', fname);
end

save('hss_underdetermined_normeq.mat', 'results', 'aspect_ratios');
fprintf('\nAll benchmarks complete.\n');


%% Visualize results
close all;

if ~exist('results', 'var')
    load('hss_underdetermined_normeq.mat', 'results', 'aspect_ratios');
end

c_matlab = [0.00 0.45 0.74];
c_normeq = [0.85 0.33 0.10];
c_build  = [0.47 0.67 0.19];

n_ar  = length(results);
n_col = ceil(n_ar / 2);
n_row = ceil(n_ar / n_col);

% -----------------------------------------------------------------------
% Figure 1: Solve time vs M
% -----------------------------------------------------------------------
figure('Name', 'Solve Time', 'Position', [50 50 1200 700]);
for ar_idx = 1:n_ar
    r     = results(ar_idx);
    M     = r.M;
    valid = ~isnan(r.matlab_lsqr_time);

    subplot(n_row, n_col, ar_idx); hold on; grid on;
    title(sprintf('AR %.1f:1', r.aspect_ratio), 'FontWeight', 'bold');
    xlabel('M (rows)'); ylabel('Avg solve time (s)');

    plot(M(valid), r.matlab_lsqr_time(valid), 'o-', 'Color', c_matlab, ...
         'LineWidth', 1.8, 'MarkerFaceColor', c_matlab, 'MarkerSize', 7, ...
         'DisplayName', 'MATLAB lsqr');
    plot(M(valid), r.hss_normeq_time(valid),  's-', 'Color', c_normeq, ...
         'LineWidth', 1.8, 'MarkerFaceColor', c_normeq, 'MarkerSize', 7, ...
         'DisplayName', 'HSS normal-eq');
    plot(M(valid), r.hss_build_time(valid),   '^--','Color', c_build, ...
         'LineWidth', 1.4, 'MarkerFaceColor', c_build,  'MarkerSize', 6, ...
         'DisplayName', 'HSS build');

    set(gca, 'XScale', 'log', 'YScale', 'log');
    legend('Location', 'northwest', 'FontSize', 8);
end
sgtitle('Solve Time vs M', 'FontSize', 13, 'FontWeight', 'bold');

% -----------------------------------------------------------------------
% Figure 2: Residual norm ||Ax-b|| vs M
% -----------------------------------------------------------------------
figure('Name', 'Residual Norm', 'Position', [80 80 1200 700]);
for ar_idx = 1:n_ar
    r     = results(ar_idx);
    M     = r.M;
    valid = ~isnan(r.matlab_residual);

    subplot(n_row, n_col, ar_idx); hold on; grid on;
    title(sprintf('AR %.1f:1', r.aspect_ratio), 'FontWeight', 'bold');
    xlabel('M (rows)'); ylabel('||Ax - b||_2');

    plot(M(valid), r.matlab_residual(valid),     'o-', 'Color', c_matlab, ...
         'LineWidth', 1.8, 'MarkerFaceColor', c_matlab, 'MarkerSize', 7, ...
         'DisplayName', 'MATLAB lsqr');
    plot(M(valid), r.hss_normeq_residual(valid), 's-', 'Color', c_normeq, ...
         'LineWidth', 1.8, 'MarkerFaceColor', c_normeq, 'MarkerSize', 7, ...
         'DisplayName', 'HSS normal-eq');

    set(gca, 'XScale', 'log', 'YScale', 'log');
    legend('Location', 'best', 'FontSize', 8);
end
sgtitle('Residual Norm ||Ax - b||_2 vs M', 'FontSize', 13, 'FontWeight', 'bold');

% -----------------------------------------------------------------------
% Figure 3: Solution error ||x_hss - x_ml|| vs M
% -----------------------------------------------------------------------
figure('Name', 'Solution Error', 'Position', [110 110 1200 700]);
for ar_idx = 1:n_ar
    r     = results(ar_idx);
    M     = r.M;
    valid = ~isnan(r.normeq_error);

    subplot(n_row, n_col, ar_idx); hold on; grid on;
    title(sprintf('AR %.1f:1', r.aspect_ratio), 'FontWeight', 'bold');
    xlabel('M (rows)'); ylabel('||x_{HSS} - x_{MATLAB}||_2');

    plot(M(valid), r.normeq_error(valid), 's-', 'Color', c_normeq, ...
         'LineWidth', 1.8, 'MarkerFaceColor', c_normeq, 'MarkerSize', 7);

    set(gca, 'XScale', 'log', 'YScale', 'log');
    grid on;
end
sgtitle('Solution Error vs MATLAB Baseline', 'FontSize', 13, 'FontWeight', 'bold');

% -----------------------------------------------------------------------
% Figure 4: ||x|| (min-norm check) vs M
% -----------------------------------------------------------------------
figure('Name', 'Solution Norm', 'Position', [140 140 1200 700]);
for ar_idx = 1:n_ar
    r     = results(ar_idx);
    M     = r.M;
    valid = ~isnan(r.matlab_xnorm);

    subplot(n_row, n_col, ar_idx); hold on; grid on;
    title(sprintf('AR %.1f:1', r.aspect_ratio), 'FontWeight', 'bold');
    xlabel('M (rows)'); ylabel('||x||_2');

    plot(M(valid), r.matlab_xnorm(valid), 'o-', 'Color', c_matlab, ...
         'LineWidth', 1.8, 'MarkerFaceColor', c_matlab, 'MarkerSize', 7, ...
         'DisplayName', 'MATLAB lsqr');
    plot(M(valid), r.hss_xnorm(valid),   's-', 'Color', c_normeq, ...
         'LineWidth', 1.8, 'MarkerFaceColor', c_normeq, 'MarkerSize', 7, ...
         'DisplayName', 'HSS normal-eq');

    set(gca, 'XScale', 'log', 'YScale', 'log');
    legend('Location', 'best', 'FontSize', 8);
end
sgtitle('Solution Norm ||x||_2 vs M  (min-norm check)', 'FontSize', 13, 'FontWeight', 'bold');

% -----------------------------------------------------------------------
% Figure 5: Speedup vs MATLAB lsqr
% -----------------------------------------------------------------------
figure('Name', 'Speedup', 'Position', [170 170 1200 700]);
for ar_idx = 1:n_ar
    r     = results(ar_idx);
    M     = r.M;
    valid = ~isnan(r.matlab_lsqr_time) & r.matlab_lsqr_time > 0 & ~isnan(r.hss_normeq_time);

    subplot(n_row, n_col, ar_idx); hold on; grid on;
    title(sprintf('AR %.1f:1', r.aspect_ratio), 'FontWeight', 'bold');
    xlabel('M (rows)'); ylabel('Speedup vs MATLAB lsqr');

    speedup = r.matlab_lsqr_time(valid) ./ r.hss_normeq_time(valid);
    plot(M(valid), speedup, 's-', 'Color', c_normeq, 'LineWidth', 1.8, ...
         'MarkerFaceColor', c_normeq, 'MarkerSize', 7);

    yline(1, '--k', 'Breakeven', 'LabelHorizontalAlignment', 'left');
    set(gca, 'XScale', 'log');
end
sgtitle('Speedup of HSS Normal-Eq vs MATLAB lsqr', 'FontSize', 13, 'FontWeight', 'bold');

fprintf('Figures rendered.\n');