%% HSS underdetermined solve benchmark:
%  Fixed aspect ratio, varying condition number and problem size
%
%  Methods compared:
%    (1) MATLAB backslash on full matrix        (if matrix fits)
%    (2) MATLAB LSQR using operator matvecs
%    (3) HSS LSQR
%    (4) HSS normal equations                  (if matrix fits)
%
%  We FIX N/M and vary:
%       - problem size
%       - kernel used to create A
%
%  Accuracy metrics:
%       ||Ax-b||_2 / ||b||_2
%       ||x-x_true||_2 / ||x_true||_2
%
%  NOTE:
%    The HSS normal-equation solve forms:
%         B = A*A'
%    explicitly before compressing into HSS.
%
% -------------------------------------------------------------------------

clear; clc;

%% ------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------

Ms             = round(logspace(3,4.2,8));   % row sizes
aspect_ratio   = 2.5;                        % FIXED aspect ratio N/M
Ns             = round(aspect_ratio * Ms);

kernel_types   = {'controlled','random','oscillatory','decay'}; 
k              = 45;
nreps          = 3;

max_full_size  = 15000;
max_bytes      = 8e9;

save_results   = true;

%% ------------------------------------------------------------------------
% Storage
% -------------------------------------------------------------------------

results = struct();

%% ========================================================================
% LOOP OVER CONDITIONING
% ========================================================================

for kidx = 1:length(kernel_types)

    kernel = kernel_types{kidx};

    fprintf('\n');
    fprintf('====================================================\n');
    %fprintf('kernel_type = %.1e\n', kernel);
    fprintf('Aspect ratio N/M = %.2f\n', aspect_ratio);
    fprintf('====================================================\n');

    % ---------------------------------------------------------------------
    % Timing
    % ---------------------------------------------------------------------

    matlab_backslash_time = nan(size(Ms));
    matlab_lsqr_time      = nan(size(Ms));
    hss_lsqr_time         = nan(size(Ms));
    hss_normeq_time       = nan(size(Ms));

    hss_build_time        = nan(size(Ms));
    hss_B_build_time      = nan(size(Ms));

    % ---------------------------------------------------------------------
    % Accuracy
    % ---------------------------------------------------------------------

    matlab_backslash_res  = nan(size(Ms));
    matlab_lsqr_res       = nan(size(Ms));
    hss_lsqr_res          = nan(size(Ms));
    hss_normeq_res        = nan(size(Ms));

    matlab_backslash_err  = nan(size(Ms));
    matlab_lsqr_err       = nan(size(Ms));
    hss_lsqr_err          = nan(size(Ms));
    hss_normeq_err        = nan(size(Ms));

    % ---------------------------------------------------------------------
    % Conditioning
    % ---------------------------------------------------------------------

    condA_vals            = nan(size(Ms));

    % ---------------------------------------------------------------------
    % Storage
    % ---------------------------------------------------------------------

    full_storage          = nan(size(Ms));
    hss_storage           = nan(size(Ms));
    gram_storage          = nan(size(Ms));
    hss_B_storage         = nan(size(Ms));

    %% ====================================================================
    % LOOP OVER PROBLEM SIZE
    % ====================================================================

    for idx = 1:length(Ms)

        M = Ms(idx);
        N = Ns(idx);

        fprintf('\n--------------------------------------------\n');
        fprintf('M = %d   N = %d   (N/M = %.2f)\n', M, N, N/M);
        fprintf('--------------------------------------------\n');

        can_form_full = (max(M,N) < max_full_size);

        %% ----------------------------------------------------------------
        % Build HSS matrix
        % -----------------------------------------------------------------

        blocksize = min(max(200, round(0.5*sqrt(M))), 2000);

        fprintf('Building HSS matrix...\n');

        [Z_handle, HZ_rect, build_time] = rectsampler(M, N, k, blocksize,'kernel_type',kernel);

        hss_build_time(idx) = build_time;

        hss_info = whos('HZ_rect');
        hss_storage(idx) = hss_info.bytes / 1e6;

        fprintf('HSS build time: %.3f s\n', build_time);

        %% ----------------------------------------------------------------
        % Full matrix if possible
        % -----------------------------------------------------------------

        if can_form_full

            fprintf('Forming full matrix...\n');

            A = full(HZ_rect);

            condA_vals(idx) = cond(A);

            Ainfo = whos('A');
            full_storage(idx) = Ainfo.bytes / 1e6;

            fprintf('cond(A) = %.3e\n', condA_vals(idx));

        else

            fprintf('Skipping full matrix formation.\n');

            full_storage(idx) = (8*M*N)/1e6;
        end

        %% ----------------------------------------------------------------
        % Create test problem with KNOWN MINIMUM-NORM SOLUTION
        % -----------------------------------------------------------------
        % 
        x_true = randn(N, 1);
        x_true = x_true / norm(x_true); % Normalize
        if can_form_full
            b = A * x_true;
        else
            b = HZ_rect *x_true;
        end
        
        % y_true = randn(M,1);
        % 
        % if can_form_full
        % 
        %     x_true = A' * y_true;
        % 
        %     % normalize to avoid scaling issues
        %     x_true = x_true / norm(x_true);
        % 
        %     b = A * x_true;
        % 
        % else
        % 
        %     % use HSS matvecs only
        %     x_true = HZ_rect' * y_true;
        % 
        %     x_true = x_true / norm(x_true);
        % 
        %     b = HZ_rect * x_true;
        % 
        % end

        %% =================================================================
        % METHOD 1:
        % MATLAB backslash (full matrix)
        % =================================================================

        if can_form_full

            fprintf('Testing MATLAB backslash...\n');

            tic;
            for r = 1:nreps
                x1 = A \ b;
            end
            matlab_backslash_time(idx) = toc / nreps;

            matlab_backslash_res(idx) = norm(A*x1 - b) / norm(b);
            matlab_backslash_err(idx) = norm(x1 - x_true) / norm(x_true);

            fprintf('  time = %.4f s\n', matlab_backslash_time(idx));
            fprintf('  rel residual = %.2e\n', matlab_backslash_res(idx));
            fprintf('  rel error    = %.2e\n', matlab_backslash_err(idx));

        else
            fprintf('Skipping MATLAB backslash.\n');
        end

        %% =================================================================
        % METHOD 2:
        % MATLAB LSQR on explicit full matrix
        % =================================================================
        
        if can_form_full
        
            fprintf('Testing MATLAB LSQR (full matrix)...\n');
        
            tic;
            for r = 1:nreps
                [x2,flag2] = lsqr(A,b,0,1000);
            end
            matlab_lsqr_time(idx) = toc / nreps;
        
            matlab_lsqr_res(idx) = norm(A*x2 - b) / norm(b);
            matlab_lsqr_err(idx) = norm(x2 - x_true) / norm(x_true);
        
            fprintf('  flag = %d\n', flag2);
            fprintf('  time = %.4f s\n', matlab_lsqr_time(idx));
            fprintf('  rel residual = %.2e\n', matlab_lsqr_res(idx));
            fprintf('  rel error    = %.2e\n', matlab_lsqr_err(idx));
        
        else
        
            fprintf('Skipping MATLAB LSQR (matrix too large).\n');
        
        end
        %% =================================================================
        % METHOD 3:
        % HSS LSQR
        % =================================================================

        fprintf('Testing HSS LSQR...\n');

        tic;
        for r = 1:nreps

            [x3,flag3] = hss_lsqr(HZ_rect,b,0,1000);
        end

        hss_lsqr_time(idx) = toc / nreps;

        if can_form_full
            hss_lsqr_res(idx) = norm(A*x3 - b) / norm(b);
        else
            hss_lsqr_res(idx) = norm(HZ_rect*x3 - b) / norm(b);
        end

        hss_lsqr_err(idx) = norm(x3 - x_true) / norm(x_true);

        fprintf('  flag = %d\n', flag3);
        fprintf('  time = %.4f s\n', hss_lsqr_time(idx));
        fprintf('  rel residual = %.2e\n', hss_lsqr_res(idx));
        fprintf('  rel error    = %.2e\n', hss_lsqr_err(idx));

        %% =================================================================
        % METHOD 4:
        % HSS normal equations
        % =================================================================

        if can_form_full

            fprintf('Testing HSS normal equations...\n');

            fprintf('  Forming B = A*A'' ...\n');

            tB = tic;
            B = A*A';
            hss_B_build_time(idx) = toc(tB);

            Binfo = whos('B');
            gram_storage(idx) = Binfo.bytes / 1e6;

            fprintf('  Compressing B into HSS...\n');

            tHB = tic;
            HB = hss(B, blocksize = blocksize);
            hss_B_build_time(idx) = ...
                hss_B_build_time(idx) + toc(tHB);

            HBinfo = whos('HB');
            hss_B_storage(idx) = HBinfo.bytes / 1e6;

            tic;
            for r = 1:nreps

                y = HB \ b;
                x4 = A' * y;

            end
            hss_normeq_time(idx) = toc / nreps;

            hss_normeq_res(idx) = norm(A*x4 - b) / norm(b);
            hss_normeq_err(idx) = norm(x4 - x_true) / norm(x_true);

            fprintf('  time = %.4f s\n', hss_normeq_time(idx));
            fprintf('  rel residual = %.2e\n', hss_normeq_res(idx));
            fprintf('  rel error    = %.2e\n', hss_normeq_err(idx));

        else

            fprintf('Skipping HSS normal equations.\n');
        end

        %% ----------------------------------------------------------------
        % Storage summary
        % -----------------------------------------------------------------

        fprintf('\nStorage summary:\n');
        fprintf('  Full matrix A     : %.1f MB\n', full_storage(idx));
        fprintf('  HSS(A)            : %.1f MB\n', hss_storage(idx));

        if can_form_full
            fprintf('  Gram matrix B     : %.1f MB\n', gram_storage(idx));
            fprintf('  HSS(B)            : %.1f MB\n', hss_B_storage(idx));
        end

        %% ----------------------------------------------------------------
        % Cleanup
        % -----------------------------------------------------------------

        clear A B HB
        clear x1 x2 x3 x4 y
        clear HZ_rect

    end

    %% ====================================================================
    % Store results
    % ====================================================================

    results(kidx).kernel                   = kernel;
    results(kidx).aspect_ratio             = aspect_ratio;

    results(kidx).M                        = Ms(:);
    results(kidx).N                        = Ns(:);

    results(kidx).condA                    = condA_vals(:);

    results(kidx).matlab_backslash_time    = matlab_backslash_time(:);
    results(kidx).matlab_lsqr_time         = matlab_lsqr_time(:);
    results(kidx).hss_lsqr_time            = hss_lsqr_time(:);
    results(kidx).hss_normeq_time          = hss_normeq_time(:);

    results(kidx).matlab_backslash_res     = matlab_backslash_res(:);
    results(kidx).matlab_lsqr_res          = matlab_lsqr_res(:);
    results(kidx).hss_lsqr_res             = hss_lsqr_res(:);
    results(kidx).hss_normeq_res           = hss_normeq_res(:);

    results(kidx).matlab_backslash_err     = matlab_backslash_err(:);
    results(kidx).matlab_lsqr_err          = matlab_lsqr_err(:);
    results(kidx).hss_lsqr_err             = hss_lsqr_err(:);
    results(kidx).hss_normeq_err           = hss_normeq_err(:);

    results(kidx).table = table( ...
        Ms(:), Ns(:), condA_vals(:), ...
        matlab_backslash_time(:), ...
        matlab_lsqr_time(:), ...
        hss_lsqr_time(:), ...
        hss_normeq_time(:), ...
        matlab_backslash_res(:), ...
        matlab_lsqr_res(:), ...
        hss_lsqr_res(:), ...
        hss_normeq_res(:), ...
        matlab_backslash_err(:), ...
        matlab_lsqr_err(:), ...
        hss_lsqr_err(:), ...
        hss_normeq_err(:), ...
        'VariableNames', { ...
            'M', 'N', 'condA', ...
            'MATLAB_backslash_s', ...
            'MATLAB_lsqr_s', ...
            'HSS_lsqr_s', ...
            'HSS_normeq_s', ...
            'MATLAB_backslash_res', ...
            'MATLAB_lsqr_res', ...
            'HSS_lsqr_res', ...
            'HSS_normeq_res', ...
            'MATLAB_backslash_err', ...
            'MATLAB_lsqr_err', ...
            'HSS_lsqr_err', ...
            'HSS_normeq_err'});

    %% Save table

    fname = sprintf('hss_kernel_type_%1.0e.txt', kernel);
    writetable(results(kidx).table, fname, 'Delimiter', '\t');

end

%% ------------------------------------------------------------------------
% Save all results
% -------------------------------------------------------------------------

if save_results
    save('hss_underdetermined_scaling.mat', ...
         'results', ...
         'kernel_types', ...
         'aspect_ratio');
end

fprintf('\nBenchmark complete.\n');
%% load
load('hss_underdetermined_scaling.mat') %or v2
%% plot_benchmark_results.m
%  Publication-quality figures for HSS solver benchmark.
%  Expects:  results(kidx)  struct as defined in the benchmark script.
%
%  Produces (per kernel):
%    Fig A – Solve time  vs N  (log-log)
%    Fig B – Relative forward error vs N  (log-log)
%    Fig C – Relative residual  vs N  (log-log)
%  And one combined summary figure with all kernels side-by-side.
%
%  Usage:
%    Run your benchmark to populate  results , then run this script.
% =========================================================================

% ── Appearance constants ─────────────────────────────────────────────────
FONT        = 'Times New Roman';   % serif font looks better in papers
FONT_SZ     = 12;
TITLE_SZ    = 13;
LABEL_SZ    = 13;
LEGEND_SZ   = 10;
LW          = 1.6;                 % line width
MS          = 6;                   % marker size
FIG_W       = 14;                  % cm  (single column ≈ 8.5, double ≈ 17)
FIG_H       = 10;

% Colour-blind-safe palette (Wong 2011, Nature Methods)
C = struct( ...
    'backslash', [0   114 178]/255, ...   % blue
    'mlsqr',     [230 159   0]/255, ...   % orange
    'hlsqr',     [0  158 115]/255, ...   % green
    'hnormeq',   [213  94   0]/255);      % vermilion

STYLES = {'-o', '--s', '-.^', ':d'};   % one per solver
SOLVER_KEYS   = {'backslash', 'mlsqr', 'hlsqr', 'hnormeq'};
SOLVER_LABELS = {'MATLAB \backslash', 'MATLAB LSQR', ...
                 'HSS LSQR', 'HSS NormEq'};

% ── Helper: apply common axis style ──────────────────────────────────────
function style_ax(ax, fs, font)
    set(ax, 'FontSize', fs, 'FontName', font, ...
            'XColor', 'black', 'YColor', 'black', ...
            'Box', 'on', 'TickDir', 'out', ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'LineWidth', 0.8, 'GridAlpha', 0.25, ...
            'GridLineStyle', '--');
    grid(ax, 'on');
end

% ── Helper: export figure ─────────────────────────────────────────────────
function export_fig_pub(fig, fname)
    exportgraphics(fig, [fname '.pdf'], ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(fig, [fname '.png'], ...
        'Resolution', 300, 'BackgroundColor', 'white');
    fprintf('  Saved  %s.{pdf,png}\n', fname);
end

% =========================================================================
%  PER-KERNEL FIGURES
% =========================================================================
nKernels = numel(results);

for kidx = 1:nKernels

    r      = results(kidx);
    M      = r.M;          % row dimension  (x-axis)
    kname  = strrep(r.kernel, '_', '\_');   % safe for TeX interpreter
    fstem  = sprintf('benchmark_%s', results(kidx).kernel);

    % ── collect data into cell arrays for looping ─────────────────────
    time_data = { r.matlab_backslash_time, r.matlab_lsqr_time, ...
                  r.hss_lsqr_time,         r.hss_normeq_time };
    res_data  = { r.matlab_backslash_res,  r.matlab_lsqr_res, ...
                  r.hss_lsqr_res,          r.hss_normeq_res };
    err_data  = { r.matlab_backslash_err,  r.matlab_lsqr_err, ...
                  r.hss_lsqr_err,          r.hss_normeq_err };

    colors = {C.backslash, C.mlsqr, C.hlsqr, C.hnormeq};

    % ==================================================================
    %  FIGURE A – Solve time
    % ==================================================================
    figA = figure('Units','centimeters','Position',[2 2 FIG_W FIG_H], ...
                  'Color','white');
    ax = axes(figA);

    hL = gobjects(4,1);
    for s = 1:4
        hL(s) = loglog(ax, M, time_data{s}, STYLES{s}, ...
            'Color', colors{s}, ...
            'LineWidth', LW, ...
            'MarkerSize', MS, ...
            'MarkerFaceColor', colors{s}, ...
            'DisplayName', SOLVER_LABELS{s});
        hold(ax,'on');
    end
    hold(ax,'off');

    style_ax(ax, FONT_SZ, FONT);
    xlabel(ax, '$M$', 'Interpreter','latex', ...
           'FontSize',LABEL_SZ, 'Color','black');
    ylabel(ax, 'Wall-clock time (s)', ...
           'FontSize',LABEL_SZ, 'FontName',FONT, 'Color','black');
    leg = legend(ax, hL, 'Location','northwest', ...
                 'FontSize',LEGEND_SZ,'FontName',FONT,'Box','on');
    leg.Title.String = 'Solver';

    add_ref_line(ax, M, time_data{1}, 2, '$\mathcal{O}(M^2)$');   % pinned to MATLAB backslash
    add_ref_line(ax, M, time_data{2}, 1, '$\mathcal{O}(M)$');     % pinned to MATLAB LSQR

    set(ax,'XTick', unique(round(logspace(log10(min(M)), ...
                                          log10(max(M)), 6))));

    export_fig_pub(figA, [fstem '_time']);

    % ==================================================================
    %  FIGURE B – Forward error
    % ==================================================================
    figB = figure('Units','centimeters','Position',[2 2 FIG_W FIG_H], ...
                  'Color','white');
    ax2 = axes(figB);

    for s = 1:4
        loglog(ax2, M, err_data{s}, STYLES{s}, ...
            'Color', colors{s}, ...
            'LineWidth', LW, ...
            'MarkerSize', MS, ...
            'MarkerFaceColor', colors{s}, ...
            'DisplayName', SOLVER_LABELS{s});
        hold(ax2,'on');
    end
    % Machine-epsilon reference line
    yline(ax2, eps, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
          'Label', '\epsilon_{mach}', 'LabelHorizontalAlignment','left', ...
          'FontSize', LEGEND_SZ, 'FontName', FONT, ...
          'Interpreter', 'tex');
    hold(ax2,'off');

    style_ax(ax2, FONT_SZ, FONT);
    xlabel(ax2, '$M$', 'Interpreter','latex', ...
           'FontSize',LABEL_SZ, 'Color','black');
    ylabel(ax2, 'Relative forward error  $\|x - x^*\| / \|x^*\|$', ...
           'Interpreter','latex','FontSize',LABEL_SZ,'Color','black');
    legend(ax2, 'Location','best','FontSize',LEGEND_SZ,'FontName',FONT,'Box','on');

    export_fig_pub(figB, [fstem '_error']);

    % ==================================================================
    %  FIGURE C – Relative residual
    % ==================================================================
    figC = figure('Units','centimeters','Position',[2 2 FIG_W FIG_H], ...
                  'Color','white');
    ax3 = axes(figC);

    for s = 1:4
        loglog(ax3, M, res_data{s}, STYLES{s}, ...
            'Color', colors{s}, ...
            'LineWidth', LW, ...
            'MarkerSize', MS, ...
            'MarkerFaceColor', colors{s}, ...
            'DisplayName', SOLVER_LABELS{s});
        hold(ax3,'on');
    end
    yline(ax3, eps, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
          'Label', '\epsilon_{mach}', 'LabelHorizontalAlignment','left', ...
          'FontSize', LEGEND_SZ, 'FontName', FONT, ...
          'Interpreter', 'tex');
    hold(ax3,'off');

    style_ax(ax3, FONT_SZ, FONT);
    xlabel(ax3, '$M$', 'Interpreter','latex', ...
           'FontSize',LABEL_SZ, 'Color','black');
    ylabel(ax3, 'Relative residual  $\|Ax - b\| / \|b\|$', ...
           'Interpreter','latex','FontSize',LABEL_SZ,'Color','black');
    legend(ax3, 'Location','best','FontSize',LEGEND_SZ,'FontName',FONT,'Box','on');

    export_fig_pub(figC, [fstem '_residual']);

    % ==================================================================
    %  FIGURE D – 3-panel summary for this kernel
    % ==================================================================
    figD = figure('Units','centimeters', ...
                  'Position',[2 2 FIG_W*1.9 FIG_H], 'Color','white');
    tiledlayout(figD, 1, 3, 'TileSpacing','compact','Padding','compact');

    panels = {time_data, err_data, res_data};
    ylabels = { 'Wall-clock time (s)', ...
                'Relative forward error', ...
                'Relative residual' };
    panel_tags = {'(a)', '(b)', '(c)'};

    for p = 1:3
        axp = nexttile;
        for s = 1:4
            loglog(axp, M, panels{p}{s}, STYLES{s}, ...
                'Color', colors{s}, 'LineWidth', LW, ...
                'MarkerSize', MS, 'MarkerFaceColor', colors{s}, ...
                'DisplayName', SOLVER_LABELS{s});
            hold(axp,'on');
        end
        if p > 1
            yline(axp, eps, ':', 'Color',[.5 .5 .5], 'LineWidth',1.2, ...
                  'HandleVisibility','off');
        end
        hold(axp,'off');
        style_ax(axp, FONT_SZ, FONT);
        xlabel(axp,'$M$','Interpreter','latex', ...
               'FontSize',LABEL_SZ,'Color','black');
        ylabel(axp, ylabels{p}, 'FontSize',LABEL_SZ,'FontName',FONT,'Color','black');
        % panel tag in lieu of title
        text(axp, 0.03, 0.97, panel_tags{p}, 'Units','normalized', ...
             'VerticalAlignment','top','FontSize',FONT_SZ,'FontName',FONT, ...
             'FontWeight','bold','Color','black');
        if p == 1
            add_ref_line(axp, M, time_data{1}, 2, '$\mathcal{O}(M^2)$');
            add_ref_line(axp, M, time_data{2}, 1, '$\mathcal{O}(M)$');
        end
    end

    % shared legend in last tile
    lgd = legend(nexttile(3), 'Location','southeast', ...
                 'FontSize',LEGEND_SZ,'FontName',FONT,'Box','on');  %#ok
    lgd.Title.String = 'Solver';

    export_fig_pub(figD, [fstem '_summary']);
    close(figD);   % keep individual figs open; close summary to save memory

end   % kernel loop


% =========================================================================
%  MULTI-KERNEL COMPARISON  (timing only – one row per kernel)
% =========================================================================
if nKernels > 1
    figE = figure('Units','centimeters', ...
                  'Position',[2 2 FIG_W*1.9 FIG_H*nKernels*0.7], ...
                  'Color','white');
    tl = tiledlayout(figE, nKernels, 3, ...
                     'TileSpacing','compact','Padding','compact'); %#ok<NASGU>

    col_titles  = {'Solve time (s)', 'Forward error', 'Residual'};

    for kidx = 1:nKernels
        r  = results(kidx);
        M  = r.M;
        kn = strrep(r.kernel,'_','\_');
        colors = {C.backslash, C.mlsqr, C.hlsqr, C.hnormeq};

        all_data = { ...
          {r.matlab_backslash_time, r.matlab_lsqr_time, r.hss_lsqr_time, r.hss_normeq_time}, ...
          {r.matlab_backslash_err,  r.matlab_lsqr_err,  r.hss_lsqr_err,  r.hss_normeq_err}, ...
          {r.matlab_backslash_res,  r.matlab_lsqr_res,  r.hss_lsqr_res,  r.hss_normeq_res} };

        for p = 1:3
            axe = nexttile;
            for s = 1:4
                loglog(axe, M, all_data{p}{s}, STYLES{s}, ...
                    'Color', colors{s}, 'LineWidth', LW, ...
                    'MarkerSize', MS-1, 'MarkerFaceColor', colors{s}, ...
                    'DisplayName', SOLVER_LABELS{s});
                hold(axe,'on');
            end
            if p > 1
                yline(axe, eps, ':', 'Color',[.5 .5 .5],'LineWidth',1.1, ...
                      'HandleVisibility','off');
            end
            hold(axe,'off');
            style_ax(axe, FONT_SZ-1, FONT);
            if kidx == nKernels
                xlabel(axe,'$M$','Interpreter','latex', ...
                       'FontSize',LABEL_SZ-1,'Color','black');
            end
            if p == 1
                ylabel(axe, kn, 'FontSize',LABEL_SZ,'FontName',FONT, ...
                       'Interpreter','tex','Color','black');
            end
            if kidx == 1
                text(axe, 0.5, 1.04, col_titles{p}, 'Units','normalized', ...
                     'HorizontalAlignment','center','FontSize',FONT_SZ, ...
                     'FontName',FONT,'Color','black','FontWeight','bold');
            end
        end
    end

    % one shared legend
    lgd2 = legend(axe,'Location','southeast','FontSize',LEGEND_SZ-1, ...
                  'FontName',FONT,'Box','on');
    lgd2.Title.String = 'Solver';

    export_fig_pub(figE,'benchmark_all_kernels');
end


% =========================================================================
%  LOCAL FUNCTION – draw a black reference slope line on a log-log axis
% =========================================================================
function add_ref_line(ax, M, y_anchor_data, slope, legend_label)
%ADD_REF_LINE  Overlay a dashed black reference line with a given log-log slope.
%   The line is pinned to pass through the geometric mean of y_anchor_data
%   at the geometric mean of M, so it tracks the nominated solver curve.
%
%   Arguments:
%     ax            – target axes
%     M             – x-data vector (row counts)
%     y_anchor_data – y-values of the solver to pin against (same length as M)
%     slope         – log-log slope (e.g. 1 for O(M), 2 for O(M^2))
%     legend_label  – string shown in the legend (LaTeX)

    x_vals   = M(:);
    y_vals   = y_anchor_data(:);
    valid    = isfinite(y_vals) & y_vals > 0;

    x_geo    = 10^mean(log10(x_vals(valid)));
    y_geo    = 10^mean(log10(y_vals(valid)));

    % y = c * x^slope  through (x_geo, y_geo)
    c        = y_geo / x_geo^slope;
    x_line   = logspace(log10(min(x_vals)), log10(max(x_vals)), 120);
    y_line   = c .* x_line .^ slope;

    hold(ax, 'on');
    plot(ax, x_line, y_line, '--k', ...
         'LineWidth',   1.1, ...
         'DisplayName', legend_label);
    hold(ax, 'off');
end

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