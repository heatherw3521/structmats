    %% Rectangular HSS underdetermined solve test via normal equations (A*A^T)
    % NOTE: Matrix-matrix HSS multiply not yet implemented.
    %       Workaround: form B = A*A' explicitly, then build HSS on B.
    %       Min-norm solve: x = A' * (A*A')^{-1} * b
    clear; clc;
    
    % -----------------------------------------------------------------------
    % Parameters
    % -----------------------------------------------------------------------
    Ms           = round(logspace(3, 4, 8));   % Row dimensions (M)
    aspect_ratios = [1.5, 2.0, 2.5, 3.0];        % N/M ratios (underdetermined: N > M)
    k            = 45;                            % HSS rank parameter
    nreps        = 3;                             % Solve repetitions for timing
    max_full_size = 12000;                        % Skip full-matrix methods above this
    
    % -----------------------------------------------------------------------
    % Storage
    % -----------------------------------------------------------------------
    results = struct();
    
    % -----------------------------------------------------------------------
    % Main loops
    % -----------------------------------------------------------------------
    for ar_idx = 1:length(aspect_ratios)
        aspect_ratio = aspect_ratios(ar_idx);
        N_sizes = round(Ms * aspect_ratio);
    
        fprintf('\n########################################\n');
        fprintf('ASPECT RATIO: %.1f:1  (N > M, underdetermined)\n', aspect_ratio);
        fprintf('########################################\n');
    
        % --- pre-allocate timing / accuracy / storage arrays ---------------
        hss_build_time      = nan(size(Ms));
        matlab_lsqr_time    = nan(size(Ms));
        hss_normeq_time     = nan(size(Ms));   % HSS solve of B=A*A', then backproject
    
        matlab_residual     = nan(size(Ms));
        hss_normeq_residual = nan(size(Ms));
    
        full_storage        = nan(size(Ms));   % MB for A (full)
        gram_storage        = nan(size(Ms));   % MB for B = A*A' (full, M×M)
        hss_storage         = nan(size(Ms));   % MB for HSS(B)
    
        % -------------------------------------------------------------------
        for idx = 1:length(Ms)
            M = Ms(idx);
            N = N_sizes(idx);
    
            fprintf('\n====================================\n');
            fprintf('Size: M = %d, N = %d  (%.1fx underdetermined)\n', M, N, N/M);
            fprintf('====================================\n');
    
            can_form_full = (M * N < max_full_size^2);
    
            % ----------------------------------------------------------------
            % Build a test matrix A  (M x N, structured/random)
            % Using a Cauchy-like construction to get low off-diagonal rank.
            % ----------------------------------------------------------------
            fprintf('Constructing test matrix A (%d x %d)...\n', M, N);
            s = (1:M)' / M;          % source points  (M x 1)
            t = (1:N)  / N;          % target points  (1 x N)
            A = 1 ./ (s + t);        % Cauchy matrix, rank-k off-diagonal blocks
    
            % ---- Form B = A * A'  (M x M square, symmetric positive semi-def)
            fprintf('Forming B = A * A''  (%d x %d)...\n', M, M);
            B = A * A';
    
            gram_info        = whos('B');
            gram_storage(idx) = gram_info.bytes / 1e6;
    
            % ---- Build HSS on B -----------------------------------------------
            blocksize = min(max(200, round(0.5*sqrt(M))), 2000);
            fprintf('Building HSS(B) with blocksize %d...\n', blocksize);
    
            t0 = tic;
            HB = hss(B, blocksize = blocksize);     % <-- adjust call to your HSS constructor
            hss_build_time(idx) = toc(t0);
    
            hss_info         = whos('HB');
            hss_storage(idx) = hss_info.bytes / 1e6;
    
            fprintf('HSS build time: %.3f s\n', hss_build_time(idx));
    
            % ---- Create test problem  b = A * x_true --------------------------
            x_true    = randn(N, 1);
            x_true    = x_true / norm(x_true);
            b         = A * x_true;               % M x 1
    
            % ================================================================
            % Method 1: MATLAB lsqr on full A  (baseline, only if A fits)
            % ================================================================
            if can_form_full
                full_info        = whos('A');
                full_storage(idx) = full_info.bytes / 1e6;
    
                fprintf('Testing MATLAB lsqr on full A...\n');
                tic;
                for r = 1:nreps
                    x_ml = lsqr(A, b, 1e-10, 2000);
                end
                matlab_lsqr_time(idx) = toc / nreps;
    
                matlab_residual(idx) = norm(A * x_ml - b);
                fprintf('Avg MATLAB lsqr: %.4f s  |  ||Ax-b|| = %.2e\n', ...
                        matlab_lsqr_time(idx), matlab_residual(idx));
            else
                full_storage(idx) = (M * N * 8) / 1e6;
                fprintf('Skipping MATLAB lsqr (A too large: %d x %d)\n', M, N);
            end
    
            % ================================================================
            % Method 2: Min-norm via HSS normal equations
            %   Solve  B * y = b  with B = A*A' (HSS)
            %   then   x = A' * y
            %
            %   This gives the min-norm least-squares solution because
            %   x = A'(AA')^{-1}b  lies in row-space of A.
            % ================================================================
            fprintf('Testing HSS normal-equations solve  (B y = b, x = A'' y)...\n');
            tic;
            for r = 1:nreps
                y_hss  = HB \ b;          % HSS solve of M x M system
                x_hss  = A' * y_hss;      % back-project: N x 1
            end
            hss_normeq_time(idx) = toc / nreps;
    
            hss_normeq_residual(idx) = norm(A * x_hss - b);
            fprintf('Avg HSS normal-eq solve: %.4f s  |  ||Ax-b|| = %.2e\n', ...
                    hss_normeq_time(idx), hss_normeq_residual(idx));
    
            % ---- Optional: check min-norm property against MATLAB solution ---
            if can_form_full
                fprintf('  ||x_hss|| = %.4e   ||x_lsqr|| = %.4e  (min-norm check)\n', ...
                        norm(x_hss), norm(x_ml));
            end
    
            % ---- Storage summary ---------------------------------------------
            fprintf('\nStorage: A (full) = %.1f MB | B = A*A'' = %.1f MB | HSS(B) = %.1f MB\n', ...
                    full_storage(idx), gram_storage(idx), hss_storage(idx));
            fprintf('Compression:  B / HSS(B) = %.1fx\n', ...
                    gram_storage(idx) / hss_storage(idx));
    
            % ---- Clean up ----------------------------------------------------
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
        results(ar_idx).full_storage        = full_storage(:);
        results(ar_idx).gram_storage        = gram_storage(:);
        results(ar_idx).hss_storage         = hss_storage(:);
    
        results(ar_idx).table = table( ...
            Ms(:), N_sizes(:), ...
            hss_build_time(:), ...
            matlab_lsqr_time(:), hss_normeq_time(:), ...
            matlab_residual(:),  hss_normeq_residual(:), ...
            full_storage(:), gram_storage(:), hss_storage(:), ...
            'VariableNames', { ...
                'M', 'N', ...
                'HSS_build_s', ...
                'MATLAB_lsqr_s', 'HSS_normeq_s', ...
                'MATLAB_residual', 'HSS_normeq_residual', ...
                'A_full_MB', 'B_gram_MB', 'HSS_B_MB'});
    
        fname = sprintf('hss_underdetermined_AR%.1f.txt', aspect_ratio);
        writetable(results(ar_idx).table, fname, 'Delimiter', '\t');
        fprintf('\nSaved results to %s\n', fname);
    end
    
    save('hss_underdetermined_normeq.mat', 'results', 'aspect_ratios');
    fprintf('\nAll benchmarks complete.\n');
    
    
    %% Visualize HSS underdetermined solve benchmark results
    clear; clc; close all;
    
    if ~exist('results', 'var')
        load('hss_underdetermined_normeq.mat', 'results', 'aspect_ratios');
    end
    
    % -----------------------------------------------------------------------
    % Style
    % -----------------------------------------------------------------------
    c_matlab = [0.00 0.45 0.74];   % blue
    c_normeq = [0.85 0.33 0.10];   % orange
    c_build  = [0.47 0.67 0.19];   % green  (build time only)
    
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
        title(sprintf('Aspect ratio %.1f:1', r.aspect_ratio), 'FontWeight', 'bold');
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
    sgtitle('Solve Time vs Problem Size', 'FontSize', 13, 'FontWeight', 'bold');
    
    % -----------------------------------------------------------------------
    % Figure 2: Residual norm ||Ax-b|| vs M
    % -----------------------------------------------------------------------
    figure('Name', 'Residual Norm', 'Position', [80 80 1200 700]);
    for ar_idx = 1:n_ar
        r     = results(ar_idx);
        M     = r.M;
        valid = ~isnan(r.matlab_residual);
    
        subplot(n_row, n_col, ar_idx); hold on; grid on;
        title(sprintf('Aspect ratio %.1f:1', r.aspect_ratio), 'FontWeight', 'bold');
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
    sgtitle('Residual Norm ||Ax - b||_2 vs Problem Size', 'FontSize', 13, 'FontWeight', 'bold');
    
    % -----------------------------------------------------------------------
    % Figure 3: Speedup of HSS normal-eq vs MATLAB lsqr
    % -----------------------------------------------------------------------
    figure('Name', 'Speedup', 'Position', [110 110 1200 700]);
    for ar_idx = 1:n_ar
        r     = results(ar_idx);
        M     = r.M;
        valid = ~isnan(r.matlab_lsqr_time) & r.matlab_lsqr_time > 0 & ~isnan(r.hss_normeq_time);
    
        subplot(n_row, n_col, ar_idx); hold on; grid on;
        title(sprintf('Aspect ratio %.1f:1', r.aspect_ratio), 'FontWeight', 'bold');
        xlabel('M (rows)'); ylabel('Speedup vs MATLAB lsqr');
    
        speedup = r.matlab_lsqr_time(valid) ./ r.hss_normeq_time(valid);
        plot(M(valid), speedup, 's-', 'Color', c_normeq, 'LineWidth', 1.8, ...
             'MarkerFaceColor', c_normeq, 'MarkerSize', 7);
    
        yline(1, '--k', 'Breakeven', 'LabelHorizontalAlignment', 'left');
        set(gca, 'XScale', 'log');
    end
    sgtitle('Speedup of HSS Normal-Eq vs MATLAB lsqr', 'FontSize', 13, 'FontWeight', 'bold');
    
    fprintf('Figures rendered.\n');