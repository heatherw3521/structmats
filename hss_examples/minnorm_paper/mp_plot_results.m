function mp_plot_results(resdir, figdir, tag)
%MP_PLOT_RESULTS  Figures of the suite from results/*.mat (re-runnable alone).
if nargin < 3, tag = ''; end
[C, MK] = mp_style();
ld = @(name) load(fullfile(resdir, [name tag '.mat']));
% ---------------- E2 stability
try
    S = ld('E2_conditioning');
    f = figure('Visible', 'off', 'Position', [100 100 900 340]);
    sets = {'graded', 'blur'}; ttl = {'(a) F2 with graded columns', '(b) F4 Gaussian blur'};
    for p = 1:2
        subplot(1, 2, p); T = S.(sets{p});
        k = [T.kappa]; [k, o] = sort(k); T = T(o);
        loglog(k, max([T.ulv_fwd], 1e-17), ['-' MK{1}], 'Color', C(1,:)); hold on
        loglog(k, [T.ne_fwd], ['-' MK{2}], 'Color', C(2,:));
        loglog(k, [T.cg_fwd], ['-' MK{3}], 'Color', C(3,:));
        loglog(k, k*eps/2, '-', 'Color', [.55 .55 .55]);
        loglog(k, k.^2*eps/2, ':', 'Color', [.55 .55 .55]);
        xlabel('\kappa(H)'); ylabel('difference from dense QR solution'); title(ttl{p}); grid on
        if p == 1
            legend('H\\b (ULV)', 'normal equations', 'CGNE', '\kappa u', '\kappa^2 u', 'Location', 'northwest');
        end
    end
    mp_figsave(f, figdir, ['E2_stability' tag]);
catch err, fprintf('E2 plot skipped: %s\n', err.message); end
% ---------------- E3/E4 scaling
try
    S = ld('E3_scaling'); R = S.rows;
    fams = {'F1','F2','F3','F4','F5'};
    f = figure('Visible', 'off', 'Position', [100 100 900 340]);
    subplot(1,2,1);
    sel = @(fam, meth) R(strcmp({R.family}, fam) & strncmp({R.method}, meth, numel(meth)));
    T = sel('F2', 'ulv'); loglog([T.n], [T.t], ['-' MK{1}], 'Color', C(1,:)); hold on
    loglog([T.n], [T.t_repeat], ['--' MK{1}], 'Color', C(1,:));
    D = sel('F2', 'dense_qr'); loglog([D.n], [D.t], ['-' MK{2}], 'Color', C(2,:));
    G = sel('F2', 'cgne'); loglog([G.n], [G.t], ['-' MK{3}], 'Color', C(3,:));
    n = [T.n]; loglog(n, T(end).t*n/n(end), '-', 'Color', [.55 .55 .55]);
    xlabel('n (columns)'); ylabel('time (s)'); grid on; title('F2, m = n/2');
    legend('H\\b, first solve', 'H\\b, repeat solve', 'dense QR', 'CGNE', 'O(n)', 'Location', 'northwest');
    subplot(1,2,2);
    for j = 1:numel(fams)
        T = sel(fams{j}, 'ulv'); semilogx([T.n], 1e6*[T.t]./[T.n], ['-' MK{j}], 'Color', C(j,:)); hold on
    end
    xlabel('n (columns)'); ylabel('first H\\b time / n  (\mus)'); grid on; title('time per column');
    legend(fams, 'Location', 'northwest');
    mp_figsave(f, figdir, ['E3_scaling' tag]);
    f = figure('Visible', 'off', 'Position', [100 100 450 340]);
    for j = 1:numel(fams)
        T = sel(fams{j}, 'ulv'); loglog([T.n], [T.err], ['-' MK{j}], 'Color', C(j,:)); hold on
    end
    xlabel('n (columns)'); ylabel('||x - x_*|| / ||x_*||'); grid on; title('manufactured min-norm error');
    legend(fams, 'Location', 'northwest');
    mp_figsave(f, figdir, ['E3_accuracy' tag]);
catch err, fprintf('E3 plot skipped: %s\n', err.message); end
% ---------------- E5/E6 sweeps
try
    S = ld('E5_sweeps');
    f = figure('Visible', 'off', 'Position', [100 100 1100 320]);
    subplot(1,3,1); bss = unique([S.rank.bs]); lg = {};
    for j = 1:numel(bss)
        R = S.rank([S.rank.bs] == bss(j));
        loglog([R.k], [R.t], ['-' MK{j}], 'Color', C(j,:)); hold on
        lg{end+1} = sprintf('blocksize %d', bss(j));
    end
    xlabel('HSS rank k'); ylabel('H\\b time (s)'); title('rank, fixed leaf size'); grid on
    legend(lg{:}, 'Location', 'northwest');
    subplot(1,3,2); semilogx([S.blocksize.bs], [S.blocksize.t], ['-' MK{2}], 'Color', C(2,:));
    xlabel('blocksize'); ylabel('H\\b time (s)'); title('blocksize'); grid on
    subplot(1,3,3); T = S.aspect;
    semilogy([T.ratio], [T.t], ['-' MK{3}], 'Color', C(3,:)); hold on
    bad = ~[T.slack_ok];
    semilogy([T(bad).ratio], [T(bad).t], 'o', 'MarkerSize', 10, 'Color', C(8,:));
    xlabel('aspect ratio m/n'); ylabel('H\\b time (s)'); title('aspect ratio (circled: Assumption 1 fails)'); grid on
    mp_figsave(f, figdir, ['E5_E6_sweeps' tag]);
catch err, fprintf('E5 plot skipped: %s\n', err.message); end
% ---------------- applications
try
    S = ld('A_applications');
    f = figure('Visible', 'off', 'Position', [100 100 900 320]);
    A1 = S.A1.gap;
    subplot(1,2,1); plot(A1.tt, A1.ftrue, '-', 'Color', [.6 .6 .6]); hold on
    plot(A1.tt, A1.fmn, '-', 'Color', C(1,:)); plot(A1.xs, A1.fx, '.', 'Color', C(2,:));
    xlabel('x'); title('A1: min-norm band-limited interpolant (gap 0.42-0.52)');
    legend('true signal', 'min-norm interpolant', 'samples'); grid on
    subplot(1,2,2); T = S.A3;
    loglog([T.n], [T.rank_T], ['-' MK{1}], 'Color', C(1,:)); hold on
    loglog([T.n], [T.rank_C], ['-' MK{2}], 'Color', C(2,:));
    xlabel('n'); ylabel('max HSS rank'); legend('Toeplitz T itself', 'Cauchy-like C = F T D^* F^*', 'Location', 'northwest');
    title('A3: ranks before / after the transform'); grid on
    mp_figsave(f, figdir, ['A_applications' tag]);
catch err, fprintf('A plot skipped: %s\n', err.message); end
% ---------------- W/T: weighted min-norm and Tikhonov
try
    S = ld('W_T');
    f = figure('Visible', 'off', 'Position', [100 100 1100 320]);
    subplot(1,3,1); T = S.W1(strcmp({S.W1.kind}, 'diag')); fams = unique({T.family});
    for j = 1:numel(fams)
        U = T(strcmp({T.family}, fams{j})); loglog([U.kappa], [U.err], ['-' MK{j}], 'Color', C(j,:)); hold on
    end
    k = [T.kappa]; loglog(sort(k), sort(k)*eps/2, '-', 'Color', [.55 .55 .55]);
    xlabel('\kappa(H L^{-1})'); ylabel('error vs dense'); legend([fams, {'\kappa u'}], 'Location', 'northwest');
    title('W1: weighted min-norm'); grid on
    subplot(1,3,2); T = S.T1([S.T1.lam_rel] > 0); shp = unique(arrayfun(@(r) sprintf('%dx%d', r.m, r.n), T, 'UniformOutput', false));
    for j = 1:numel(shp)
        U = T(strcmp(arrayfun(@(r) sprintf('%dx%d', r.m, r.n), T, 'UniformOutput', false), shp{j}));
        loglog([U.lam_rel], [U.err], ['-' MK{j}], 'Color', C(j,:)); hold on
    end
    xlabel('\lambda / ||H||'); ylabel('error vs SVD filter'); legend(shp, 'Location', 'northeast'); title('T1: Tikhonov'); grid on
    subplot(1,3,3); T = S.T2(strcmp({S.T2.shape}, 'wide')); U = S.T2(strcmp({S.T2.shape}, 'tall'));
    loglog([T.n], [T.t_minnorm], ['-' MK{1}], 'Color', C(1,:)); hold on
    loglog([T.n], [T.t_tikhonov], ['-' MK{2}], 'Color', C(2,:));
    loglog([U.n], [U.t_tikhonov], ['-' MK{3}], 'Color', C(3,:));
    xlabel('n'); ylabel('time (s)'); legend('min-norm H\\b (wide)', 'Tikhonov (wide)', 'Tikhonov (tall, m = 2n)', 'Location', 'northwest');
    title('T2: cost'); grid on
    mp_figsave(f, figdir, ['W_T' tag]);
catch err, fprintf('W/T plot skipped: %s\n', err.message); end
try
    S = ld('apps2');
    f = figure('Visible', 'off', 'Position', [100 100 1100 320]);
    subplot(1,3,1); semilogy(S.BP.irls(:,1), S.BP.irls(:,3), ['-' MK{1}], 'Color', C(1,:)); hold on
    semilogy(S.BP.dr(:,1), S.BP.dr(:,3), ['-' MK{2}], 'Color', C(2,:));
    xlabel('iteration'); ylabel('||x_k - x_0|| / ||x_0||'); legend('IRLS (weighted min-norm)', 'Douglas-Rachford (projection)'); title('BP: sparse spikes'); grid on
    subplot(1,3,2); R = S.DB.rows; loglog(R(:,1), R(:,2), ['-' MK{1}], 'Color', C(1,:)); hold on
    xline_ = S.DB.lam_dp; loglog([xline_ xline_], [min(R(:,2)) max(R(:,2))], '-', 'Color', [.55 .55 .55]);
    xlabel('\lambda'); ylabel('||x_\lambda - x_0|| / ||x_0||'); title('DB: Tikhonov deblurring (line: discrepancy)'); grid on
    subplot(1,3,3); plot(S.MRI.x0, '-', 'Color', [.6 .6 .6]); hold on
    plot(S.MRI.x_mn, '-', 'Color', C(2,:)); plot(S.MRI.x_w, '-', 'Color', C(1,:));
    legend('phantom', 'min-norm', 'k-space weighted'); title('MRI: non-Cartesian 1-D'); grid on
    mp_figsave(f, figdir, ['apps2' tag]);
catch err, fprintf('apps2 plot skipped: %s\n', err.message); end
end
