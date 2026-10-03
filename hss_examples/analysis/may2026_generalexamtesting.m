%% ========================================================================
% HSS RECTANGULAR SCALING BENCHMARK
%
% Goal:
%   Separate:
%       (A) Construction + matvec scaling
%       (B) Solve scaling
%
% Construction/matvec tests scale to:
%       M,N ~ O(1e5)
%
% Solve tests remain moderate-sized because:
%       - dense MATLAB methods become infeasible
%       - HSS normal equations require forming A*A'
%
% Problem:
%       Rectangular underdetermined matrices
%       N = aspect_ratio * M
%
% Matrix family:
%       Structured rectangular kernels from rectsampler()
%
% Outputs:
%       - .mat results
%       - publication-quality figures
%
% ========================================================================

clear; clc; close all;

%% ========================================================================
% PARAMETERS
% ========================================================================

aspect_ratio = 2.5;

kernel_types = {'controlled','random','oscillatory','decay'};

k = 45;

nreps = 3;

save_results = true;

% ------------------------------------------------------------------------
% Large-scale matvec/construction tests
% ------------------------------------------------------------------------

Ms_large = round(logspace(3,4.5,9));   % up to 100000
Ns_large = round(aspect_ratio * Ms_large);

% ------------------------------------------------------------------------
% Solve tests (smaller)
% ------------------------------------------------------------------------

Ms_solve = round(logspace(3,4.1,7));
Ns_solve = round(aspect_ratio * Ms_solve);

max_full_size = 15000;

%% ========================================================================
% STORAGE
% ========================================================================

results = struct();

%% ========================================================================
% MAIN LOOP OVER KERNEL TYPES
% ========================================================================

for kidx = 1:length(kernel_types)

    kernel = kernel_types{kidx};

    fprintf('\n=================================================\n');
    fprintf('Kernel type: %s\n', kernel);
    fprintf('Aspect ratio: %.2f\n', aspect_ratio);
    fprintf('=================================================\n');

    %% ====================================================================
    % PART I:
    % LARGE-SCALE CONSTRUCTION + MATVEC
    % ====================================================================

    fprintf('\n');
    fprintf('=================================================\n');
    fprintf('PART I: CONSTRUCTION + MATVEC SCALING\n');
    fprintf('=================================================\n');

    build_time        = nan(size(Ms_large));
    hss_matvec_time   = nan(size(Ms_large));
    dense_matvec_time = nan(size(Ms_large));

    matvec_error      = nan(size(Ms_large));

    hss_storage       = nan(size(Ms_large));
    dense_storage     = nan(size(Ms_large));

    for idx = 1:length(Ms_large)

        M = Ms_large(idx);
        N = Ns_large(idx);

        fprintf('\n[MATVEC TEST] M = %d, N = %d\n', M, N);

        blocksize = min(max(500, round(2*sqrt(M))), 4000);

        %% ----------------------------------------------------------------
        % Build HSS
        % -----------------------------------------------------------------

        [~, H, tbuild] = rectsampler( ...
            M, N, k, blocksize, ...
            'kernel_type', kernel);

        build_time(idx) = tbuild;

        hinfo = whos('H');
        hss_storage(idx) = hinfo.bytes / 1e6;

        fprintf('  HSS build time: %.3f s\n', tbuild);

        %% ----------------------------------------------------------------
        % Matvec timing
        % -----------------------------------------------------------------

        x = randn(N,1);

        tic;
        for r = 1:nreps
            y_hss = H*x;
        end
        hss_matvec_time(idx) = toc / nreps;

        fprintf('  HSS matvec: %.4e s\n', hss_matvec_time(idx));

        %% ----------------------------------------------------------------
        % Dense matvec ONLY for moderate sizes
        % -----------------------------------------------------------------

        can_form_full = (M <= 15000);

        if can_form_full

            fprintf('  Forming dense matrix...\n');

            A = full(H);

            ainfo = whos('A');
            dense_storage(idx) = ainfo.bytes / 1e6;

            tic;
            for r = 1:nreps
                y_dense = A*x;
            end
            dense_matvec_time(idx) = toc / nreps;

            matvec_error(idx) = ...
                norm(y_hss - y_dense) / norm(y_dense);

            fprintf('  Dense matvec: %.4e s\n', dense_matvec_time(idx));
            fprintf('  Relative error: %.2e\n', matvec_error(idx));

            clear A y_dense

        else

            dense_storage(idx) = (8*M*N)/1e6;

            fprintf('  Dense matrix skipped.\n');
        end

        clear H y_hss x

    end
    results(kidx).kernel = kernel;

    % ------------------------------------------------------------
    % Large-scale matvec
    % ------------------------------------------------------------

    results(kidx).matvec.M              = Ms_large(:);
    results(kidx).matvec.N              = Ns_large(:);

    results(kidx).matvec.build_time     = build_time(:);
    results(kidx).matvec.hss_mv_time    = hss_matvec_time(:);
    results(kidx).matvec.dense_mv_time  = dense_matvec_time(:);

    results(kidx).matvec.mv_error       = matvec_error(:);

    results(kidx).matvec.hss_storage    = hss_storage(:);
    results(kidx).matvec.dense_storage  = dense_storage(:);

    continue
    
    %% ====================================================================
    % PART II:
    % SOLVE TESTS
    % ====================================================================

    fprintf('\n');
    fprintf('=================================================\n');
    fprintf('PART II: SOLVE BENCHMARKS\n');
    fprintf('=================================================\n');

    matlab_lsqminnorm_time = nan(size(Ms_solve));
    matlab_lsqr_time      = nan(size(Ms_solve));
    hss_lsqr_time         = nan(size(Ms_solve));
    hss_normeq_time       = nan(size(Ms_solve));

    matlab_lsqminnorm_res  = nan(size(Ms_solve));
    matlab_lsqr_res       = nan(size(Ms_solve));
    hss_lsqr_res          = nan(size(Ms_solve));
    hss_normeq_res        = nan(size(Ms_solve));

    matlab_lsqminnorm_err  = nan(size(Ms_solve));
    matlab_lsqr_err       = nan(size(Ms_solve));
    hss_lsqr_err          = nan(size(Ms_solve));
    hss_normeq_err        = nan(size(Ms_solve));

    condA_vals            = nan(size(Ms_solve));

    for idx = 1:length(Ms_solve)

        M = Ms_solve(idx);
        N = Ns_solve(idx);

        fprintf('\n[SOLVE TEST] M = %d, N = %d\n', M, N);

        blocksize = min(max(200, round(0.5*sqrt(M))), 2000);

        %% ----------------------------------------------------------------
        % Build HSS
        % -----------------------------------------------------------------

        [~, H, ~] = rectsampler( ...
            M, N, k, blocksize, ...
            'kernel_type', kernel);

        can_form_full = (max(M,N) < max_full_size);

        %% ----------------------------------------------------------------
        % Dense matrix
        % -----------------------------------------------------------------

        if can_form_full

            A = full(H);

            condA_vals(idx) = cond(A);

        end

        %% ----------------------------------------------------------------
        % Construct MINIMUM-NORM test problem
        %
        % x_true = A' y
        % b = A x_true
        % ensures x_true in row-space(A)
        % -----------------------------------------------------------------

        y_true = randn(M,1);

        if can_form_full

            x_true = A' * y_true;

            x_true = x_true / norm(x_true);

            b = A * x_true;

        else

            x_true = H' * y_true;

            x_true = x_true / norm(x_true);

            b = H * x_true;

        end

        %% =================================================================
        % METHOD 1:
        % MATLAB BACKSLASH
        % =================================================================

        if can_form_full

            fprintf('  MATLAB backslash...\n');

            tic;
            for r = 1:nreps
                x1 = A\b;
            end
            matlab_lsqminnorm_time(idx) = toc / nreps;

            matlab_lsqminnorm_res(idx) = ...
                norm(A*x1 - b) / norm(b);

            matlab_lsqminnorm_err(idx) = ...
                norm(x1 - x_true) / norm(x_true);

        end

        %% =================================================================
        % METHOD 2:
        % MATLAB LSQR (FULL MATRIX)
        % =================================================================

        if can_form_full

            fprintf('  MATLAB LSQR...\n');

            tic;
            for r = 1:nreps
                [x2,~] = lsqr(A,b,1e-12,2000);
            end
            matlab_lsqr_time(idx) = toc / nreps;

            matlab_lsqr_res(idx) = ...
                norm(A*x2 - b) / norm(b);

            matlab_lsqr_err(idx) = ...
                norm(x2 - x_true) / norm(x_true);

        end

        %% =================================================================
        % METHOD 3:
        % HSS LSQR
        % =================================================================

        fprintf('  HSS LSQR...\n');

        tic;
        for r = 1:nreps
            [x3,~] = hss_lsqr(H,b,1e-12,2000);
        end
        hss_lsqr_time(idx) = toc / nreps;

        if can_form_full
            hss_lsqr_res(idx) = norm(A*x3 - b) / norm(b);
        else
            hss_lsqr_res(idx) = norm(H*x3 - b) / norm(b);
        end

        hss_lsqr_err(idx) = ...
            norm(x3 - x_true) / norm(x_true);

        %% =================================================================
        % METHOD 4:
        % HSS NORMAL EQUATIONS
        % =================================================================

        if can_form_full

            fprintf('  HSS normal equations...\n');

            B = A*A';

            HB = hss(B, blocksize = blocksize);

            tic;
            for r = 1:nreps

                y = HB\b;
                x4 = A'*y;

            end
            hss_normeq_time(idx) = toc / nreps;

            hss_normeq_res(idx) = ...
                norm(A*x4 - b) / norm(b);

            hss_normeq_err(idx) = ...
                norm(x4 - x_true) / norm(x_true);

            clear B HB y x4

        end

        clear H A

    end

    %% ====================================================================
    % STORE RESULTS
    % ====================================================================

    results(kidx).kernel = kernel;

    % ------------------------------------------------------------
    % Large-scale matvec
    % ------------------------------------------------------------

    results(kidx).matvec.M              = Ms_large(:);
    results(kidx).matvec.N              = Ns_large(:);

    results(kidx).matvec.build_time     = build_time(:);
    results(kidx).matvec.hss_mv_time    = hss_matvec_time(:);
    results(kidx).matvec.dense_mv_time  = dense_matvec_time(:);

    results(kidx).matvec.mv_error       = matvec_error(:);

    results(kidx).matvec.hss_storage    = hss_storage(:);
    results(kidx).matvec.dense_storage  = dense_storage(:);

    % ------------------------------------------------------------
    % Solve
    % ------------------------------------------------------------

    % results(kidx).solve.M = Ms_solve(:);
    % results(kidx).solve.N = Ns_solve(:);
    % 
    % results(kidx).solve.condA = condA_vals(:);
    % 
    % results(kidx).solve.matlab_backslash_time = matlab_backslash_time(:);
    % results(kidx).solve.matlab_lsqr_time      = matlab_lsqr_time(:);
    % results(kidx).solve.hss_lsqr_time         = hss_lsqr_time(:);
    % results(kidx).solve.hss_normeq_time       = hss_normeq_time(:);
    % 
    % results(kidx).solve.matlab_backslash_res = matlab_backslash_res(:);
    % results(kidx).solve.matlab_lsqr_res      = matlab_lsqr_res(:);
    % results(kidx).solve.hss_lsqr_res         = hss_lsqr_res(:);
    % results(kidx).solve.hss_normeq_res       = hss_normeq_res(:);
    % 
    % results(kidx).solve.matlab_backslash_err = matlab_backslash_err(:);
    % results(kidx).solve.matlab_lsqr_err      = matlab_lsqr_err(:);
    % results(kidx).solve.hss_lsqr_err         = hss_lsqr_err(:);
    % results(kidx).solve.hss_normeq_err       = hss_normeq_err(:);

end

%% ========================================================================
% SAVE RESULTS
% ========================================================================

if save_results

    save('hss_rectangular_scaling_results.mat', ...
        'results', ...
        'kernel_types', ...
        'aspect_ratio');

end

fprintf('\nAll experiments complete.\n');

%% ========================================================================
% PUBLISHABLE FIGURES
% ========================================================================

set(groot,'defaultAxesFontSize',14);
set(groot,'defaultAxesLineWidth',1.2);
set(groot,'defaultLineLineWidth',2);
set(groot,'defaultFigureColor','w');

cols = lines(4);

%% Load Results
load("hss_rectangular_scaling_results.mat")

%% plot_matvec_results.m
%  Publication-quality figures for HSS matvec / storage / build-time benchmark.
%  Expects:  results  struct loaded from hss_rectangular_scaling_results.mat
%
%  Figures produced:
%    Fig 1  – HSS construction time vs M          (log-log)
%    Fig 2  – Storage: HSS vs dense vs M          (log-log)
%    Fig 3  – Matvec time: HSS vs dense vs M      (log-log, ref lines)
%    Fig 4  – Matvec relative error vs M          (semi-log y)
% =========================================================================

% ── Appearance constants (identical to plot_benchmark_results.m) ──────────
FONT      = 'Times New Roman';
FONT_SZ   = 12;
LABEL_SZ  = 13;
LEGEND_SZ = 10;
LW        = 1.6;
MS        = 6;
FIG_W     = 14;    % cm
FIG_H     = 10;

% Wong (2011) colour-blind-safe palette – one colour per kernel
KCOLS = [ ...
      0  114  178;   % blue
    230  159    0;   % orange
      0  158  115;   % green
    213   94    0;   % vermilion
     86  180  233;   % sky blue
    204  121  167;   % reddish purple
    240  228   66] / 255;  % yellow

% ── Helpers ───────────────────────────────────────────────────────────────
function style_ax(ax, fs, font)
    set(ax, 'FontSize', fs, 'FontName', font, ...
            'XColor', 'black', 'YColor', 'black', ...
            'Box', 'on', 'TickDir', 'out', ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'LineWidth', 0.8, 'GridAlpha', 0.25, ...
            'GridLineStyle', '--');
    grid(ax, 'on');
end

function export_fig_pub(fig, fname)
    exportgraphics(fig, [fname '.pdf'], ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(fig, [fname '.png'], ...
        'Resolution', 300, 'BackgroundColor', 'white');
    fprintf('  Saved  %s.{pdf,png}\n', fname);
end

function add_ref_line(ax, M, y_anchor_data, slope, legend_label)
    x_vals = M(:);
    y_vals = y_anchor_data(:);
    valid  = isfinite(y_vals) & y_vals > 0;
    x_geo  = 10^mean(log10(x_vals(valid)));
    y_geo  = 10^mean(log10(y_vals(valid)));
    c      = y_geo / x_geo^slope;
    x_line = logspace(log10(min(x_vals)), log10(max(x_vals)), 120);
    y_line = c .* x_line .^ slope;
    hold(ax, 'on');
    plot(ax, x_line, y_line, '--k', ...
         'LineWidth', 1.1, 'DisplayName', legend_label);
    hold(ax, 'off');
end

% ── Load data ─────────────────────────────────────────────────────────────
load('hss_rectangular_scaling_results.mat');
nK = numel(results);

% =========================================================================
%  FIGURE 1 – HSS construction time
% =========================================================================
fig1 = figure('Units','centimeters','Position',[2 2 FIG_W FIG_H],'Color','white');
ax1  = axes(fig1);
hold(ax1, 'on');

ci = 1;
for kidx = 2:nK
    r = results(kidx);
    loglog(ax1, r.matvec.M, r.matvec.build_time, '-o', ...
        'Color',           KCOLS(ci,:), ...
        'MarkerFaceColor', KCOLS(ci,:), ...
        'MarkerSize',      MS, ...
        'LineWidth',       LW, ...
        'DisplayName',     r.kernel);

    ci = ci + 1;
end

hold(ax1, 'off');
style_ax(ax1, FONT_SZ, FONT);
xlabel(ax1, '$M$', 'Interpreter','latex','FontSize',LABEL_SZ,'Color','black');
ylabel(ax1, 'Construction time (s)', 'FontName',FONT,'FontSize',LABEL_SZ,'Color','black');
legend(ax1, 'Location','northwest','FontSize',LEGEND_SZ,'FontName',FONT,'Box','on');

%export_fig_pub(fig1, 'may2026hss_buildtime');

% =========================================================================
%  FIGURE 2 – Storage: HSS vs dense
% =========================================================================
fig2 = figure('Units','centimeters','Position',[2 2 FIG_W FIG_H],'Color','white');
ax2  = axes(fig2);


ci = 1;
for kidx = 4:nK
    r     = results(kidx);
    valid = ~isnan(r.matvec.dense_storage);

    loglog(ax2, r.matvec.M, r.matvec.hss_storage, '-o', ...
        'Color',           KCOLS(ci,:), ...
        'MarkerFaceColor', KCOLS(ci,:), ...
        'MarkerSize',      MS, ...
        'LineWidth',       LW, ...
        'DisplayName',     [r.kernel ' (HSS)']);
    hold(ax2, 'on');
    loglog(ax2, r.matvec.M(valid), r.matvec.dense_storage(valid), '--s', ...
        'Color',           KCOLS(ci,:), ...
        'MarkerFaceColor', 'white', ...
        'MarkerSize',      MS, ...
        'LineWidth',       LW, ...
        'DisplayName',     [r.kernel ' (dense)']);
    
    ci = ci + 1;
end

hold(ax2, 'off');
style_ax(ax2, FONT_SZ, FONT);
xlabel(ax2, '$M$', 'Interpreter','latex','FontSize',LABEL_SZ,'Color','black');
ylabel(ax2, 'Storage ($10^6$ bits)', 'Interpreter','latex','FontSize',LABEL_SZ,'Color','black');
legend(ax2, 'Location','northwest','FontSize',LEGEND_SZ,'FontName',FONT,'Box','on');

%export_fig_pub(fig2, 'may2026hss_storage');

% =========================================================================
%  FIGURE 3 – Matvec timing: HSS vs dense  (+ reference lines)
% =========================================================================
fig3 = figure('Units','centimeters','Position',[2 2 FIG_W FIG_H],'Color','white');
ax3  = axes(fig3);


ci = 1;
for kidx = 2:nK
    r     = results(kidx);
    valid = ~isnan(r.matvec.dense_mv_time);

    loglog(ax3, r.matvec.M, r.matvec.hss_mv_time, '-o', ...
        'Color',           KCOLS(ci,:), ...
        'MarkerFaceColor', KCOLS(ci,:), ...
        'MarkerSize',      MS, ...
        'LineWidth',       LW, ...
        'DisplayName',     [r.kernel ' (HSS)']);
    hold(ax3, 'on');
    loglog(ax3, r.matvec.M(valid), r.matvec.dense_mv_time(valid), '--s', ...
        'Color',           KCOLS(ci,:), ...
        'MarkerFaceColor', 'white', ...
        'MarkerSize',      MS, ...
        'LineWidth',       LW, ...
        'DisplayName',     [r.kernel ' (dense)']);
    ci = ci + 1;
end

% Reference lines – O(M) pinned to HSS, O(M^2) pinned to dense
r_ref   = results(2);
valid2  = ~isnan(r_ref.matvec.dense_mv_time);
add_ref_line(ax3, r_ref.matvec.M,          r_ref.matvec.hss_mv_time,           1, '$\mathcal{O}(M)$');
add_ref_line(ax3, r_ref.matvec.M(valid2),  r_ref.matvec.dense_mv_time(valid2), 2, '$\mathcal{O}(M^2)$');

hold(ax3, 'off');
style_ax(ax3, FONT_SZ, FONT);
xlabel(ax3, '$M$', 'Interpreter','latex','FontSize',LABEL_SZ,'Color','black');
ylabel(ax3, 'Time (s)', 'FontName',FONT,'FontSize',LABEL_SZ,'Color','black');
legend(ax3, 'Location','northwest','FontSize',LEGEND_SZ,'FontName',FONT,'Box','on', ...
       'Interpreter','latex');

%export_fig_pub(fig3, 'may2026matvec_scaling');

% =========================================================================
%  FIGURE 4 – Matvec relative error
% =========================================================================
fig4 = figure('Units','centimeters','Position',[2 2 FIG_W FIG_H],'Color','white');
ax4  = axes(fig4);


ci = 1;
for kidx = 1:nK
    r     = results(kidx);
    valid = ~isnan(r.matvec.mv_error);

    semilogy(ax4, r.matvec.M(valid), r.matvec.mv_error(valid), '-o', ...
        'Color',           KCOLS(ci,:), ...
        'MarkerFaceColor', KCOLS(ci,:), ...
        'MarkerSize',      MS, ...
        'LineWidth',       LW, ...
        'DisplayName',     r.kernel);
    hold(ax4, 'on');
    ci = ci + 1;
end

% Machine-epsilon floor
yline(ax4, eps, ':', 'Color',[0.5 0.5 0.5], 'LineWidth',1.1, ...
      'Label','\epsilon_{mach}', 'LabelHorizontalAlignment','left', ...
      'FontSize',LEGEND_SZ, 'FontName',FONT, 'Interpreter','tex', ...
      'HandleVisibility','off');

hold(ax4, 'off');
style_ax(ax4, FONT_SZ, FONT);
xlabel(ax4, '$M$', 'Interpreter','latex','FontSize',LABEL_SZ,'Color','black');
ylabel(ax4, 'Relative matvec error', 'FontName',FONT,'FontSize',LABEL_SZ,'Color','black');
legend(ax4, 'Location','southwest','FontSize',LEGEND_SZ,'FontName',FONT,'Box','on');

%export_fig_pub(fig4, 'may2026matvec_accuracy');


%% ========================================================================
% FIGURE 1:
% HSS Build Time scaling
% ========================================================================

figure('Position',[100 100 900 650]);

hold on;

for kidx = 2:length(results)

    r = results(kidx);

    loglog( ...
        r.matvec.M, ...
        r.matvec.build_time, ...
        '-o', ...
        'Color', cols(kidx,:), ...
        'MarkerFaceColor', cols(kidx,:), ...
        'DisplayName', r.kernel);

end

grid on;
box on;

xlabel('$M$','Interpreter','latex');
ylabel('Construction time (s)','Interpreter','latex');
%title('HSS Construction Scaling','Interpreter','latex');

legend('Location','northwest');

exportgraphics(gcf,'may2026hss_buildtime.pdf','ContentType','vector');

%% ========================================================================
% FIGURE 1.5:
% Storage scaling
% ========================================================================

figure('Position',[100 100 900 650]);



for kidx = 4:length(results)

    r = results(kidx);

    loglog( ...
        r.matvec.M, ...
        r.matvec.hss_storage, ...
        '-o', ...
        'Color', cols(kidx,:), ...
        'MarkerFaceColor', cols(kidx,:), ...
        'DisplayName', "HSS");
    hold on;
    valid = ~isnan(r.matvec.dense_mv_time);
    loglog( ...
        r.matvec.M(valid), ...
        r.matvec.dense_storage(valid), ...
        '--s', ...
        'Color', cols(kidx,:), ...
        'MarkerFaceColor', 'w', ...
        'DisplayName', ['Dense']);
end

grid on;
box on;

xlabel('$M$','Interpreter','latex');
ylabel('Storage ($1e6*$ bits)','Interpreter','latex');
%title('HSS Construction Scaling','Interpreter','latex');

legend('Location','northwest');

exportgraphics(gcf,'may2026hss_storage.pdf','ContentType','vector');

%% ========================================================================
% FIGURE 2: Matvec scaling (HSS vs dense)
% ========================================================================

figure('Position',[100 100 900 650]); 
ax = axes;
set(ax, 'XScale', 'log', 'YScale', 'log', ...
        'Box', 'off', ...                          % no top/right border
        'XGrid', 'off', 'YGrid', 'on', ...         % horizontal gridlines only
        'GridLineStyle', '--', ...
        'GridColor', [0.7 0.7 0.7], ...            % light grey dashes
        'GridAlpha', 1.0, ...
        'TickDir', 'out', ...
        'FontSize', 13, ...
        'TickLabelInterpreter', 'latex');

for kidx = 2:length(results)

    r = results(kidx);

    % -------------------------
    % HSS matvec
    % -------------------------
    loglog( ...
        r.matvec.M, ...
        r.matvec.hss_mv_time, ...
        '-o', ...
        'Color', cols(kidx,:), ...
        'MarkerFaceColor', cols(kidx,:), ...
        'MarkerSize',      6, ...
        'LineWidth',       1.5, ...
        'DisplayName', [r.kernel ' (HSS)']);
    hold on;
    % -------------------------
    % Dense matvec (only valid entries)
    % -------------------------
    valid = ~isnan(r.matvec.dense_mv_time);

    loglog( ...
        r.matvec.M(valid), ...
        r.matvec.dense_mv_time(valid), ...
        '--s', ...
        'Color', cols(kidx,:), ...
        'MarkerFaceColor', 'w', ...
        'MarkerSize',      6, ...
        'LineWidth',       1.5, ...
        'DisplayName', [r.kernel ' (dense)']);

end

% -------------------------
% Reference lines
% -------------------------
% Gather all M values and HSS/dense times across all results for anchoring
all_M          = results(2).matvec.M;
all_hss_time   = results(2).matvec.hss_mv_time;
valid2         = ~isnan(results(2).matvec.dense_mv_time);
all_dense_time = results(2).matvec.dense_mv_time(valid2);
all_M_dense    = results(2).matvec.M(valid2);

% Anchor O(n) ref line to the midpoint of the HSS data
M_ref_hss   = all_M;
mid_hss     = round(numel(M_ref_hss)/2);
c_linear    = all_hss_time(mid_hss) / all_M(mid_hss);        % t = c * n
ref_linear  = c_linear * M_ref_hss;

% Anchor O(n^2) ref line to the midpoint of the dense data
M_ref_dense = all_M_dense;
mid_dense   = round(numel(M_ref_dense)/2);
c_quad      = all_dense_time(mid_dense) / all_M_dense(mid_dense).^2; % t = c * n^2
ref_quad    = c_quad * M_ref_dense.^2;

% Plot reference lines
loglog(M_ref_hss,   ref_linear, 'k-',  'LineWidth', 1.2, ...
    'DisplayName', 'O(n) reference');
loglog(M_ref_dense, ref_quad,   'k--', 'LineWidth', 1.2, ...
    'DisplayName', 'O(n^2) reference');


grid on;

xlabel('$M$',           'Interpreter', 'latex', 'FontSize', 15);
ylabel('Time (s)',      'Interpreter', 'latex', 'FontSize', 15);

leg = legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12); 

exportgraphics(gcf,'may2026matvec_scaling.pdf','ContentType','vector');

%% ========================================================================
% FIGURE 3: Matvec accuracy
% ========================================================================

figure('Position',[100 100 900 650]); hold on;

for kidx = 1:length(results)

    r = results(kidx);

    valid = ~isnan(r.matvec.mv_error);

    semilogy( ...
        r.matvec.M(valid), ...
        r.matvec.mv_error(valid), ...
        '-o', ...
        'Color', cols(kidx,:), ...
        'MarkerFaceColor', cols(kidx,:), ...
        'DisplayName', r.kernel);

end

grid on; box on;

xlabel('$M$','Interpreter','latex');
ylabel('Relative matvec error','Interpreter','latex');

legend('Location','southwest');

exportgraphics(gcf,'may2026matvec_accuracy.pdf','ContentType','vector');

%% ========================================================================
% FIGURE 4:
% Solve timing
% ========================================================================
% 
% figure('Position',[100 100 1000 700]);
% 
% tiledlayout(2,2,'Padding','compact');
% 
% names = { ...
%     'matlab\_backslash\_time', ...
%     'matlab\_lsqr\_time', ...
%     'hss\_lsqr\_time', ...
%     'hss\_normeq\_time'};
% 
% titles = { ...
%     'MATLAB Backslash', ...
%     'MATLAB LSQR', ...
%     'HSS LSQR', ...
%     'HSS Normal Equations'};
% 
% for j = 1:4
% 
%     nexttile;
%     hold on;
% 
%     for kidx = 1:length(results)
% 
%         r = results(kidx);
% 
%         field = strrep(names{j},'\_','_');
% 
%         vals = r.solve.(field);
% 
%         valid = ~isnan(vals);
% 
%         loglog( ...
%             r.solve.M(valid), ...
%             vals(valid), ...
%             '-o', ...
%             'Color', cols(kidx,:), ...
%             'MarkerFaceColor', cols(kidx,:), ...
%             'DisplayName', r.kernel);
% 
%     end
% 
%     grid on;
%     box on;
% 
%     xlabel('$M$','Interpreter','latex');
%     ylabel('Solve time (s)','Interpreter','latex');
% 
%     title(titles{j},'Interpreter','latex');
% 
% end
% 
% legend(kernel_types,'Location','bestoutside');
% 
% exportgraphics(gcf,'may2026solve_scaling.pdf','ContentType','vector');


%% LSQR Testing

aspect_ratio = 2.5;
conds = [1e3, 1e5, 1e7, 1e9, 1e11];
conds = [10];
k = 45;

nreps = 2;

save_results = true;

% ------------------------------------------------------------------------
% Solve tests (smaller)
% ------------------------------------------------------------------------

Ms_solve = round(logspace(3,3.9,6));
Ns_solve = round(aspect_ratio * Ms_solve);

max_full_size = 15000;

results = struct();


for cind = 1:length(conds)

    condnum = conds(cind);

    fprintf('\n=================================================\n');
    fprintf('Condition number: %s\n', condnum);
    fprintf('=================================================\n');


    matlab_lsqminnorm_time = nan(size(Ms_solve));
    matlab_lsqr_time      = nan(size(Ms_solve));
    hss_lsqr_time         = nan(size(Ms_solve));
    hss_normeq_time       = nan(size(Ms_solve));

    matlab_lsqminnorm_res  = nan(size(Ms_solve));
    matlab_lsqr_res       = nan(size(Ms_solve));
    hss_lsqr_res          = nan(size(Ms_solve));
    hss_normeq_res        = nan(size(Ms_solve));

    matlab_lsqminnorm_err  = nan(size(Ms_solve));
    matlab_lsqr_err       = nan(size(Ms_solve));
    hss_lsqr_err          = nan(size(Ms_solve));
    hss_normeq_err        = nan(size(Ms_solve));


    for idx = 1:length(Ms_solve)

        M = Ms_solve(idx);
        N = Ns_solve(idx);

        fprintf('\n[SOLVE TEST] M = %d, N = %d\n', M, N);

        blocksize = min(max(200, round(0.5*sqrt(M))), 2000);
        depth = floor(log(max(M,N)/blocksize)/log(2)); % Amount of levels

        % Build HSS
        % [A,info]=hss_rect_test(M,N,k,condnum,'depth',20,'verify',true);
        % 
        % condnum = info.cond_actual
        % 
        % HA = hss(A, blocksize=blocksize, sizeA=[M, N]);

        A = hankel_rect(M,N);

        HA = hss(A,blocksize = blocksize,sizeA=[M,N]);
        norm(full(HA)-A)/norm(A)

        condnum = cond(A)

        % [A_handle, HA, varargout] = rectsampler(M, N, k, blocksize);
        % A = A_handle(1:M,1:N);
        % condnum = cond(A)

        % Construct true minnorm solution
        y_true = randn(M,1);
        x_true = A'*y_true;
        x_true = x_true/norm(x_true);
        b = A * x_true;
        
        
        % backslash
        fprintf('  MATLAB lsqminnorm...\n');

        tic;
        for r = 1:nreps
            x1 = lsqminnorm(A,b);
        end
        matlab_lsqminnorm_time(idx) = toc / nreps;

        matlab_lsqminnorm_res(idx) = ...
            norm(A*x1 - b) / norm(b);

        matlab_lsqminnorm_err(idx) = ...
            norm(x1 - x_true) / norm(x_true)
        
        % matlab lsqr
        fprintf('  MATLAB LSQR...\n');

        tic;
        for r = 1:nreps
            [x2,~] = lsqr(A,b,0,2000);
        end
        matlab_lsqr_time(idx) = toc / nreps;

        matlab_lsqr_res(idx) = ...
            norm(A*x2 - b) / norm(b);

        matlab_lsqr_err(idx) = ...
            norm(x2 - x_true) / norm(x_true)

        % hss lsqr
        fprintf('  HSS LSQR...\n');
        fun = @(x, mode) lsqr_wrapper(x, mode, HA);
        tic;
        for r = 1:nreps
            %[x3,~] = hss_lsqr(HA,b,1e-12,2000);
            [x3,~] = lsqr(fun, b, 0, 2000);
        end
        hss_lsqr_time(idx) = toc / nreps;

        hss_lsqr_res(idx) = norm(A*x3 - b) / norm(b);

        hss_lsqr_err(idx) = ...
            norm(x3 - x_true) / norm(x_true)

        % normal equations
        fprintf('  HSS normal equations...\n');

        B = A*A';

        HB = hss(B, blocksize = blocksize);

        tic;
        for r = 1:nreps

            y = HB\b;
            x4 = HA.'*y;

        end
        hss_normeq_time(idx) = toc / nreps;

        hss_normeq_res(idx) = ...
            norm(A*x4 - b) / norm(b);

        hss_normeq_err(idx) = ...
            norm(x4 - x_true) / norm(x_true)

        clear B HB y

        clear HA A

    end

    % ====================================================================
    % STORE RESULTS
    % ====================================================================

    results(cind).condnum = condnum;

    results(cind).solve.M = Ms_solve(:);
    results(cind).solve.N = Ns_solve(:);

    results(cind).solve.matlab_backslash_time = matlab_lsqminnorm_time(:);
    results(cind).solve.matlab_lsqr_time      = matlab_lsqr_time(:);
    results(cind).solve.hss_lsqr_time         = hss_lsqr_time(:);
    results(cind).solve.hss_normeq_time       = hss_normeq_time(:);

    results(cind).solve.matlab_backslash_res = matlab_lsqminnorm_res(:);
    results(cind).solve.matlab_lsqr_res      = matlab_lsqr_res(:);
    results(cind).solve.hss_lsqr_res         = hss_lsqr_res(:);
    results(cind).solve.hss_normeq_res       = hss_normeq_res(:);

    results(cind).solve.matlab_backslash_err = matlab_lsqminnorm_err(:);
    results(cind).solve.matlab_lsqr_err      = matlab_lsqr_err(:);
    results(cind).solve.hss_lsqr_err         = hss_lsqr_err(:);
    results(cind).solve.hss_normeq_err       = hss_normeq_err(:);

end

save('hss_udsolvetiming_may2026_singlecond.mat', ...
        'results', ...
        'conds');

%% Plotting LSQR testing
set(groot,'defaultAxesFontSize',14);
set(groot,'defaultAxesLineWidth',1.2);
set(groot,'defaultLineLineWidth',2);
set(groot,'defaultFigureColor','w');

cols = lines(length(results));

figure('Position',[100 100 1000 700]);

tiledlayout(2,2,'Padding','compact');

names = { ...
    'matlab\_backslash\_time', ...
    'matlab\_lsqr\_time', ...
    'hss\_lsqr\_time', ...
    'hss\_normeq\_time'};

titles = { ...
    'MATLAB LSQMinnorm', ...
    'MATLAB LSQR', ...
    'HSS LSQR', ...
    'HSS Normal Equations'};

for j = 1:4

    nexttile;
    hold on;

    for cidx = 1:length(results)

        r = results(cidx);

        field = strrep(names{j},'\_','_');

        vals = r.solve.(field);

        valid = ~isnan(vals);

        loglog( ...
            r.solve.M(valid), ...
            vals(valid), ...
            '-o', ...
            'Color', cols(cidx,:), ...
            'MarkerFaceColor', cols(cidx,:));

    end

    grid on;
    box on;

    xlabel('$M$','Interpreter','latex');
    ylabel('Solve time (s)','Interpreter','latex');

    title(titles{j},'Interpreter','latex');

end

legend(kernel_types,'Location','bestoutside');

exportgraphics(gcf,'may2026_udsolvers.pdf','ContentType','vector');



%% FIGURE: Least-squares solve errors — stacked by condition number
% ========================================================================

% --- method definitions --------------------------------------------------
%   field    : suffix in results(cind).solve.*_err
%   label    : legend entry (LaTeX)
%   marker   : line + marker spec
%   color    : RGB

meth = struct( ...
    'field',  { 'matlab_backslash_err', ...
                'matlab_lsqr_err',      ...
                'hss_lsqr_err',         ...
                'hss_normeq_err'        }, ...
    'label',  { 'MATLAB \texttt{lsqminnorm}', ...
                'MATLAB \texttt{lsqr}',        ...
                'HSS-LSQR',                    ...
                'HSS normal eq.'               }, ...
    'marker', { '-o', '-s', '-^', '-d'  }, ...
    'color',  { [0.00 0.45 0.70], ...   % blue
                [0.85 0.33 0.10], ...   % red-orange
                [0.47 0.67 0.19], ...   % green
                [0.49 0.18 0.56]  }  ...% purple
);

% --- layout --------------------------------------------------------------
nconds = numel(results);
ncols  = min(nconds, 3);               % wrap to at most 3 columns
nrows  = ceil(nconds / ncols);

fig = figure('Position', [100 100 310*ncols 380*nrows], 'Color', 'w');
tl  = tiledlayout(nrows, ncols, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- one tile per condition number ---------------------------------------
for cind = 1:nconds
    r  = results(cind);
    M  = r.solve.M;

    ax = nexttile;
    set(ax, ...
        'XScale',            'log',    ...
        'YScale',            'log',    ...
        'Box',               'off',    ...   % left + bottom axes only
        'XGrid',             'off',    ...
        'YGrid',             'on',     ...   % horizontal guide lines only
        'GridLineStyle',     '--',     ...
        'GridColor',         [0.7 0.7 0.7], ...
        'GridAlpha',         1.0,      ...
        'TickDir',           'out',    ...
        'FontSize',          11,       ...
        'TickLabelInterpreter', 'latex');
    hold(ax, 'on');

    % --- plot each method ------------------------------------------------
    for mi = 1:numel(meth)
        err   = r.solve.(meth(mi).field);
        valid = ~isnan(err) & err > 0;
        if ~any(valid), continue; end

        loglog(ax, M(valid), err(valid), meth(mi).marker, ...
            'Color',           meth(mi).color, ...
            'MarkerFaceColor', meth(mi).color, ...
            'MarkerSize',      5,              ...
            'LineWidth',       1.5,            ...
            'DisplayName',     meth(mi).label);
    end

    % --- machine epsilon reference ---------------------------------------
    yline(ax, eps, '--', ...
        'Color',             [0.5 0.5 0.5], ...
        'LineWidth',         0.8,           ...
        'HandleVisibility',  'off');

    text(ax, max(M), eps * 1.8, '$\varepsilon_{\rm mach}$', ...
        'Interpreter',       'latex',         ...
        'FontSize',          9,               ...
        'Color',             [0.5 0.5 0.5],  ...
        'HorizontalAlignment','right');

    % --- panel title: kappa = 10^exp ------------------------------------
    lg10k = log10(r.condnum);
    if abs(lg10k - round(lg10k)) < 0.05
        ttl = sprintf('$\\kappa = 10^{%d}$', round(lg10k));
    else
        ttl = sprintf('$\\kappa \\approx 10^{%.1f}$', lg10k);
    end
    title(ax, ttl, 'Interpreter', 'latex', 'FontSize', 12);

    xlabel(ax, '$M$',            'Interpreter', 'latex', 'FontSize', 13);
    ylabel(ax, 'Relative error', 'Interpreter', 'latex', 'FontSize', 13);
end

% --- shared legend in a south tile --------------------------------------
%   Collect line handles from the first tile (all methods appear there)
first_ax = nexttile(1);
h = findobj(first_ax, 'Type', 'Line', '-not', 'HandleVisibility', 'off');
h = flipud(h);      % restore draw order (findobj returns reversed)

leg = legend(first_ax, h, {meth.label}, ...
    'Interpreter',  'latex',        ...
    'FontSize',     11,             ...
    'Orientation',  'horizontal',   ...
    'NumColumns',   numel(meth));
leg.Box         = 'off';
leg.Layout.Tile = 'south';

exportgraphics(fig, 'may2026solve_errors.pdf', 'ContentType', 'vector');
%% 

[A,info]=hss_rect_test(1000,2500,45,1000,'depth',10,'verify',true);
info.cond_actual

H = hss(A);

norm(full(H)-A)/norm(A)

%% Helper functions


function [Z_handle, HZ, varargout] = rectsampler(M, N, k, blocksize, varargin)
    % Create structured low-rank rectangular matrix
    % Optional: kernel_type = 'random' (default), 'oscillatory', 'decay'
    
    p = inputParser;
    addParameter(p, 'kernel_type', 'random', @ischar);
    addParameter(p, 'freq', 10, @isnumeric);  % for oscillatory
    addParameter(p, 'decay_rate', 0.1, @isnumeric);  % for decay
    addParameter(p, 'alpha', 1e-8, @isnumeric);
    parse(p, varargin{:});
    
    kernel_type = p.Results.kernel_type;
    
    % Create geometry (like Hankel example)
    t_row = linspace(0, 1, M+1); t_row(end) = [];
    t_col = linspace(0, 1, N+1); t_col(end) = [];
    
    % Store geometry for block function
    geom.t_row = t_row(:);
    geom.t_col = t_col(:);
    geom.k = k;
    geom.kernel_type = kernel_type;
    geom.freq = p.Results.freq;
    geom.decay_rate = p.Results.decay_rate;
    geom.alpha = p.Results.alpha;
    
    % Create function handle
    Z_handle = @(I, J) rect_kernel_block(I, J, geom);
    
    % Build HSS from function handle
    tic;
    HZ = hss(Z_handle, blocksize=blocksize, sizeA=[M, N]);
    timing = toc;
    
    % Apply diagonal modification
    % [HZ,dH] = diagmodify(HZ);
    % 
    % Z_handle = @(I,J) rect_kernel_block(I,J,geom) + dH(I,J);
    Z_handle = @(I,J) rect_kernel_block(I,J,geom);

    % Return timing if requested
    if nargout == 3
        varargout{1} = timing;
    elseif nargout == 4
        varargout{1} = timing;
        varargout{2} = dH;
    end
end

function B = rect_kernel_block(I, J, geom)
    % Compute block of low-rank structured matrix
    
    ti = geom.t_row(I);
    tj = geom.t_col(J);
    k = geom.k;
    
    switch geom.kernel_type
        case 'controlled'
            rng(12345);
            U_full = randn(length(geom.t_row), k);
            V_full = randn(length(geom.t_col), k);
            
            % orthonormalize
            [U_full,~] = qr(U_full,0);
            [V_full,~] = qr(V_full,0);
            
            % prescribed singular value decay
            sigma = geom.alpha.^linspace(0,1,k);
            
            S = diag(sigma);
            
            B = U_full(I,:) * S * V_full(J,:)';
        case 'random'
            % Random low-rank (fastest, simplest)
            % Generate deterministically based on indices for consistency
            rng(12345);  % Fixed seed for reproducibility
            U_full = randn(length(geom.t_row), k) / sqrt(k);
            V_full = randn(length(geom.t_col), k) / sqrt(k);
            B = U_full(I, :) * V_full(J, :)';
            
        case 'oscillatory'
            % Oscillatory kernel (like discretized integral operator)
            % B(i,j) = sum_{l=1}^k sin(freq*l*ti) * cos(freq*l*tj)
            freq = geom.freq;
            U = zeros(length(I), k);
            V = zeros(length(J), k);
            for l = 1:k
                U(:, l) = sin(freq * l * ti) / sqrt(k);
                V(:, l) = cos(freq * l * tj) / sqrt(k);
            end
            B = U * V';
            
        case 'decay'
            % Exponentially decaying kernel
            % Mimics distance-based decay
            [Ti, Tj] = meshgrid(tj, ti);
            dist = abs(Ti - Tj);
            
            % Low-rank approximation via truncated SVD of decay kernel
            % For small blocks, just compute directly
            if length(I)*length(J) < 1000
                B = exp(-geom.decay_rate * dist);
            else
                % For larger blocks, use low-rank structure
                U = zeros(length(I), k);
                V = zeros(length(J), k);
                for l = 1:k
                    lambda = geom.decay_rate * l;
                    U(:, l) = exp(-lambda * ti) / sqrt(k);
                    V(:, l) = exp(-lambda * tj) / sqrt(k);
                end
                B = U * V';
            end
            
        case 'polynomial'
            % Polynomial kernel: (1 + ti*tj')^k approximation
            U = zeros(length(I), k);
            V = zeros(length(J), k);
            for l = 1:k
                U(:, l) = ti.^(l-1) / sqrt(factorial(l-1));
                V(:, l) = tj.^(l-1) / sqrt(factorial(l-1));
            end
            B = U * V';
            
        otherwise
            error('Unknown kernel type: %s', geom.kernel_type);
    end
end

function [H,dH] = diagmodify(H)
if H.isleaf
    dH = rand(size(H.D));
    H.D = H.D+dH;
    % [Q,~] = qr(randn(size(H.D)));
    % dg = 20.^(-0:size(H.D,1)-1);
    % blk = Q*diag(dg)*Q';
    % newblk = zeros(size(H.D));
    % newblk(1:size(blk,1),1:size(blk,2)) = blk;
    % H.D = H.D + newblk;
else
    [H.A11,dH1] = diagmodify(H.A11);
    [H.A22,dH2] = diagmodify(H.A22);
    dH = blkdiag(dH1,dH2);
end
end

function [x,flag] = hss_lsqr(HZ_rect,b,tol,maxit)

fun = @(x, mode) lsqr_wrapper(x, mode, HZ_rect);

[x,flag] = lsqr(fun, b, tol, maxit);

end

function y = lsqr_wrapper(x, mode, H)
    if strcmp(mode, 'notransp')
        y = H * x; % Ax
    else
        y = H.' * x; % A'x
    end
end

function [A, geom] = hankel_rect(M, N, varargin)
%HANKEL_RECT  Rectangular Helmholtz BIE kernel matrix with HSS structure.
%
%   [A, GEOM] = HANKEL_RECT(M, N) returns an M-by-N complex matrix
%
%       A(i,j) = (i/4) * H_0^(1)(ka * r_ij) * w_j
%
%   where r_ij = |z_tgt(i) - z_src(j)| is the target-source distance and
%   w_j is the j-th arc-length quadrature weight on the source circle.
%
%   By default the M target points lie on a unit circle and the N source
%   points lie on a concentric circle of radius r_src = 0.5.  Because the
%   two point sets are geometrically separated the kernel H_0^(1) is
%   smooth everywhere, so every off-diagonal block in the HSS tree is
%   numerically low-rank.
%
%   The effective HSS rank grows roughly as  O(ka * arc_length / level)
%   per tree level.  To target a specific HSS rank k choose ka so that
%       ka * (2*pi * min(r_tgt, r_src) / 2^depth) ~ k
%
%   Optional name-value arguments
%   ------------------------------
%   'ka'        (default 10)    Helmholtz wavenumber (dimensionless)
%   'r_tgt'     (default 1.0)   radius of target circle
%   'r_src'     (default 0.5)   radius of source circle
%                               Separation = |r_tgt - r_src|.
%                               Larger separation -> lower HSS rank.
%                               Set equal to r_tgt only when M == N and
%                               you want the original square formulation.
%   'self_int'  (default true)  When M == N and r_tgt == r_src, replace
%                               the singular diagonal with 1 (matching
%                               the original Hankel_mat convention).
%
%   Returns
%   -------
%   A     M-by-N complex matrix
%   GEOM  struct: z_tgt (M x 1), z_src (N x 1), w_src (N x 1),
%                 r (M x N pairwise distances), ka
%
%   Geometry choices and HSS quality
%   ---------------------------------
%   Separated circles (r_tgt ~= r_src)  -- recommended
%       No singularity anywhere.  All blocks, including the first-level
%       off-diagonal, are smooth -> best HSS compressibility.
%
%   Same circle, M ~= N  (r_tgt == r_src)
%       Near-diagonal entries are large (log singularity) but for M ~= N
%       the "diagonal" does not align with any single index -> no
%       automatic fix.  Use only if your HSS library handles this.
%
%   Same circle, M == N  (square, 'self_int' true)
%       Recovers the original Hankel_mat behaviour exactly.
%
%   Examples
%   --------
%   % Default: 512 targets on unit circle, 256 sources on r=0.5 circle
%   [A, g] = hankel_rect(512, 256);
%
%   % Higher wavenumber -> higher HSS rank
%   A = hankel_rect(1024, 512, 'ka', 50);
%
%   % Wider separation -> lower HSS rank (more compressible)
%   A = hankel_rect(1024, 512, 'ka', 10, 'r_src', 0.1);
%
%   % Recover original square matrix
%   A = hankel_rect(N, N, 'r_src', 1.0, 'self_int', true);

% ---- parse arguments ----------------------------------------------------
p = inputParser;
addRequired(p,  'M',        @(x) isscalar(x) && x >= 1);
addRequired(p,  'N',        @(x) isscalar(x) && x >= 1);
addParameter(p, 'ka',       10,   @(x) isscalar(x) && x > 0);
addParameter(p, 'r_tgt',    1.0,  @(x) isscalar(x) && x > 0);
addParameter(p, 'r_src',    0.5,  @(x) isscalar(x) && x > 0);
addParameter(p, 'self_int', true, @(x) islogical(x) || x==0 || x==1);
parse(p, M, N, varargin{:});
opt = p.Results;

% ---- target points (M x 1) on circle of radius r_tgt -------------------
t_tgt  = linspace(0, 2*pi, M + 1);  t_tgt(end) = [];
z_tgt  = opt.r_tgt * exp(1i * t_tgt(:));

% ---- source points (N x 1) on circle of radius r_src -------------------
t_src  = linspace(0, 2*pi, N + 1);  t_src(end) = [];
z_src  = opt.r_src * exp(1i * t_src(:));

% Arc-length quadrature weights: w_j = (2pi/N) * |dz/dt|_j
%   For a circle of radius r: |dz/dt| = r, so w_j = 2*pi*r_src/N
zp_src = 1i * opt.r_src * exp(1i * t_src(:));   % tangent vector (N x 1)
w_src  = (2*pi / N) * abs(zp_src);              % N x 1

% ---- pairwise distances: M x N -----------------------------------------
%   d(i,j) = z_tgt(i) - z_src(j)
d = bsxfun(@minus, z_tgt, z_src.');   % M x N
r = abs(d);                            % M x N, all positive when separated

% ---- kernel matrix ------------------------------------------------------
%   A(i,j) = (i/4) * H_0^(1)(ka * r_ij) * w_j
%
%   w_src (N x 1) -> w_src.' (1 x N) broadcasts across rows so that
%   each column j (source j) is multiplied by its own quadrature weight.
%   This is the correct BIE convention for a rectangular matrix.
A = (0.25i * besselh(0, opt.ka * r)) .* w_src.';   % M x N

% ---- diagonal regularisation (square, same circle only) ----------------
same_circle = (opt.r_tgt == opt.r_src);
if opt.self_int && (M == N) && same_circle
    idx    = sub2ind([M, N], 1:M, 1:M);
    A(idx) = 1;
end

% ---- output geometry ---------------------------------------------------
if nargout > 1
    geom.z_tgt = z_tgt;          % M x 1  target points
    geom.z_src = z_src;          % N x 1  source points
    geom.w_src = w_src;          % N x 1  quadrature weights
    geom.r     = r;              % M x N  pairwise distances
    geom.ka    = opt.ka;
    geom.separation = abs(opt.r_tgt - opt.r_src);
end

end % hankel_rect

function [A, info] = hss_rect_test(m, n, k, kappa, varargin)
%HSS_RECT_TEST  Rectangular HSS test matrix with controlled condition number.
%
%   [A, INFO] = HSS_RECT_TEST(M, N, K, KAPPA) returns an M-by-N matrix
%   whose off-diagonal blocks (in a balanced binary partition of rows and
%   columns) have rank K, and whose 2-norm condition number is approximately
%   KAPPA.
%
%   Parameters
%   ----------
%   M, N     : matrix dimensions
%   K        : HSS rank  (rank of every off-diagonal block)
%   KAPPA    : target 2-norm condition number (>= 1)
%
%   Optional name-value arguments
%   ------------------------------
%   'depth'      (default 1)
%       Number of recursive HSS levels.  At depth 0 the whole matrix is
%       a single leaf.  At depth 1 there are four blocks, two diagonal and
%       two off-diagonal.  Deeper trees make the structure richer but the
%       condition number bound compounds (use 'verify' to check).
%
%   'off_scale'  (default 0.1, must be in (0,1))
%       Each off-diagonal block has 2-norm  off_scale / kappa.
%       By Weyl's inequality (see Theory below) this guarantees
%           cond(A)  <=  kappa / (1 - off_scale)
%       at depth 1.  Decrease off_scale for tighter control; the default
%       gives <= 11% overshoot.
%
%   'sv_dist'    ('geometric' | 'linear',  default 'geometric')
%       How the r = min(m,n) leaf-block singular values are distributed
%       across [1/kappa, 1]:
%           geometric  --  geometrically spaced  (1, kappa^{-1/(r-1)}, ..., 1/kappa)
%           linear     --  linearly spaced        (1, ..., 1/kappa)
%       Geometric spacing is more realistic for kernel matrices.
%
%   'verify'     (default false)
%       If true, compute svd(A) and store the exact condition number in
%       INFO.cond_actual.  Avoid for large M, N.
%
%   'seed'       (default [])
%       If non-empty, call rng(seed) before any random draws for
%       reproducibility.
%
%   Returns
%   -------
%   A            M-by-N double matrix
%   INFO         struct with fields:
%       .kappa_target   requested condition number
%       .cond_bound     guaranteed upper bound  kappa / (1 - off_scale)
%                       (tight at depth 1; informative at deeper levels)
%       .cond_actual    svd-based condition number, or NaN if verify=false
%       .sigma_range    [1/kappa, 1]  -- leaf singular value interval
%       .off_scale      value used
%       .depth          value used
%       .hss_rank       K
%
%   Theory
%   ------
%   At depth 1 the matrix has the 2-by-2 block structure
%
%       A = [ A11   A12 ]     A_diag = [ A11  0  ]   A_off = [ 0    A12 ]
%           [ A21   A22 ]              [  0  A22 ]            [ A21   0  ]
%
%   Leaf blocks A11, A22 are constructed with singular values in
%   [1/kappa, 1], so sigma_min(A_diag) = 1/kappa, sigma_max(A_diag) = 1.
%
%   Off-diagonal blocks are built as
%       A12 = s * U1 * B12 * V2',   s = off_scale/kappa
%   where U1, V2 are random thin orthonormal matrices and B12 is a random
%   k-by-k matrix normalised to unit 2-norm, giving ||A12||_2 = s.
%
%   A key identity: ||A_off||_2 = max(||A12||_2, ||A21||_2) = s
%   (block off-diagonal with two diagonal zero blocks; its squared
%   singular values are the union of those of A12*A12' and A21*A21').
%
%   Weyl's inequality then gives
%       sigma_min(A) >= sigma_min(A_diag) - ||A_off||_2
%                    = 1/kappa - off_scale/kappa
%                    = (1 - off_scale) / kappa
%       sigma_max(A) <= 1 + off_scale/kappa  ~  1   (large kappa)
%   hence
%       cond(A) <= kappa / (1 - off_scale).
%
%   Example
%   -------
%   % 2-level HSS matrix, 1024 x 512, rank 8, condition number ~1e4
%   [A, info] = hss_rect_test(1024, 512, 8, 1e4, ...
%       'depth', 2, 'verify', true);
%   fprintf('target kappa = %.2e,  actual = %.2e\n', ...
%       info.kappa_target, info.cond_actual);

% ---- parse arguments ----------------------------------------------------
p = inputParser;
addRequired(p, 'm',     @(x) isscalar(x) && x >= 1);
addRequired(p, 'n',     @(x) isscalar(x) && x >= 1);
addRequired(p, 'k',     @(x) isscalar(x) && x >= 1);
addRequired(p, 'kappa', @(x) isscalar(x) && x >= 1);
addParameter(p, 'depth',     1,           @(x) isscalar(x) && x >= 0);
addParameter(p, 'off_scale', 0.1,         @(x) isscalar(x) && x > 0 && x < 1);
addParameter(p, 'sv_dist',   'geometric', @(x) ismember(x,{'geometric','linear'}));
addParameter(p, 'verify',    false,       @(x) islogical(x) || x==0 || x==1);
addParameter(p, 'seed',      [],          @(x) isempty(x) || isscalar(x));
parse(p, m, n, k, kappa, varargin{:});
opt = p.Results;

if ~isempty(opt.seed)
    rng(opt.seed);
end

% ---- build matrix -------------------------------------------------------
A = build_level(m, n, k, kappa, opt.off_scale, opt.sv_dist, opt.depth);

% ---- output info --------------------------------------------------------
info.kappa_target = kappa;
info.cond_bound   = kappa / (1 - opt.off_scale);
info.sigma_range  = [1/kappa, 1];
info.off_scale    = opt.off_scale;
info.depth        = opt.depth;
info.hss_rank     = k;
info.cond_actual  = NaN;

if opt.verify
    sv = svd(A);
    info.cond_actual = sv(1) / sv(end);
end

end % hss_rect_test

function A = build_level(m, n, k, kappa, off_scale, sv_dist, depth)
%BUILD_LEVEL  Recursively assemble one level of the HSS tree.

m1 = floor(m/2);  m2 = m - m1;
n1 = floor(n/2);  n2 = n - n1;

% Base case: leaf block (also catches m or n too small to split further)
if depth == 0 || m1 < k || n1 < k || m2 < k || n2 < k
    A = make_leaf(m, n, kappa, sv_dist);
    return;
end

% ---- diagonal children (recursive) -------------------------------------
A11 = build_level(m1, n1, k, kappa, off_scale, sv_dist, depth - 1);
A22 = build_level(m2, n2, k, kappa, off_scale, sv_dist, depth - 1);

% ---- off-diagonal blocks: rank k, 2-norm = off_scale/kappa -------------
%
%   Constructed as  U * B * V'  where
%     U (mi x k), V (nj x k)  are random thin orthonormal factors
%     B (k x k)                is random with unit 2-norm
%   => 2-norm of the block = off_scale/kappa  (exactly)
%
off_norm = off_scale / kappa;

[U1, ~] = qr(randn(m1, k), 'econ');
[V2, ~] = qr(randn(n2, k), 'econ');
B12 = randn(k, k);
B12 = B12 / norm(B12, 2);               % normalise to unit 2-norm
A12 = off_norm*(U1 * B12 * V2');      % m1 x n2, rank k, ||A12||_2 = off_norm

[U2, ~] = qr(randn(m2, k), 'econ');
[V1, ~] = qr(randn(n1, k), 'econ');
B21 = randn(k, k);
B21 = B21 / norm(B21, 2);
A21 = off_norm*(U2 * B21 * V1');      % m2 x n1, rank k, ||A21||_2 = off_norm

% ---- assemble ----------------------------------------------------------
A = [A11,  A12;
     A21,  A22];

end % build_level

function D = make_leaf(m, n, kappa, sv_dist)
%MAKE_LEAF  Full-rank m-by-n block with singular values in [1/kappa, 1].

r = min(m, n);

if r == 1 || kappa == 1
    sv = ones(1, r);
elseif strcmp(sv_dist, 'geometric')
    % sv(1) = 1,  sv(r) = 1/kappa,  geometric progression
    sv = (1/kappa) .^ ( (0 : r-1) / (r - 1) );
else
    sv = linspace(1, 1/kappa, r);
end

[U, ~] = qr(randn(m, r), 'econ');
[V, ~] = qr(randn(n, r), 'econ');
D = U * diag(sv) * V';

end % make_leaf