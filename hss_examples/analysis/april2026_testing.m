%%
ns = 2.^(8:12);   % matrix size
ls = ns/8;        % leaf size     
k = 50;           % k
fullt = [];
componentt = [];
for i=1:5
    n = ns(i);
    l = ls(i);
    A = randn(n);
    Gamma = rand(n,k);
    tic
    S = A * Gamma;
    fullt(i) = toc;

    tic
    for j=1:8
        Sj = A((j-1)*l+1:j*l,:)*Gamma;
        S = [S; Sj];
    end
    componentt(i) = toc;
end




%% testing gaussian mult with full A vs pieces of A
ns = 2.^(8:14);   % matrix sizes
k = 50;           
ntrials = 3;

fullt = zeros(size(ns));
componentt = zeros(size(ns));

for i = 1:length(ns)
    n = ns(i);
    l = n/8;

    A = randn(n);
    Gamma = randn(n,k);

    % --- optionally remove diagonal ---
    A = A - diag(diag(A));

    % --- FULL MULTIPLY ---
    t = 0;
    for ttrial = 1:ntrials
        tic
        S_full = A * Gamma;
        t = t + toc;
    end
    fullt(i) = t / ntrials;

    % --- COMPONENT (BLOCK ROWS) ---
    t = 0;
    for ttrial = 1:ntrials
        S_comp = zeros(n,k);   % preallocate

        tic
        for j = 1:8
            I = (j-1)*l + 1 : j*l;
            S_comp(I,:) = A(I,:) * Gamma;
        end
        t = t + toc;
    end
    componentt(i) = t / ntrials;

    fprintf('n = %d done\n', n);
end

loglog(ns, fullt, 'o-', 'LineWidth', 2); hold on;
loglog(ns, componentt, 's-', 'LineWidth', 2);

% reference O(n^2 k)
c = fullt(1) / (ns(1)^2);
loglog(ns, c*ns.^2, '--');

legend('Full A*Gamma', 'Block rows', 'O(n^2)', 'Location','northwest');
xlabel('n'); ylabel('time');
grid on;


%% testing gaussian vs sparsestack


function S = sparsesign_slow(d,n,zeta)
    cols = kron((1:n)',ones(zeta,1)); % zeta nonzeros per column
    vals = 2*randi(2,n*zeta,1) - 3; % uniform random +/-1 values
    rows = zeros(n*zeta,1);
    for i = 1:n
       rows((i-1)*zeta+1:i*zeta) = randsample(d,zeta);
    end
    S = sparse(rows, cols, vals / sqrt(zeta), d, n);
end

% -----------------------------
% parameters
% -----------------------------
ns = round(logspace(3, 4.2, 10));
l_frac = 0.05;
zeta = 8;
ntrials = 5;

t_gauss = zeros(length(ns),1);
t_sparse_matlab = zeros(length(ns),1);
t_sparse_mex = zeros(length(ns),1);

% -----------------------------
% main loop
% -----------------------------
for ii = 1:length(ns)

    n = ns(ii);
    l = max(10, round(l_frac * n));

    fprintf('n = %d, l = %d\n', n, l);

    A = randn(n,n);

    tg = zeros(ntrials,1);
    tsm = zeros(ntrials,1);
    tmx = zeros(ntrials,1);

    for t = 1:ntrials

        % -------------------------
        % Gaussian
        % -------------------------
        tic;
        Omega = randn(n,l) / sqrt(l);
        Y = A * Omega;
        tg(t) = toc;

        % -------------------------
        % MATLAB SparseStack
        % -------------------------
        tic;
        S = sparsesign_slow(n,l,zeta);
        Y = A * S;
        tsm(t) = toc;

        % -------------------------
        % MEX SparseStack
        % -------------------------
        tic;
        S = sparsesign(n,l,zeta);
        Y = A * S;
        tmx(t) = toc;

    end

    t_gauss(ii) = mean(tg);
    t_sparse_matlab(ii) = mean(tsm);
    t_sparse_mex(ii) = mean(tmx);

end

% -----------------------------
% plot
% -----------------------------
figure;
loglog(ns, t_gauss, '-o', 'LineWidth', 2); hold on;
loglog(ns, t_sparse_matlab, '-s', 'LineWidth', 2);
loglog(ns, t_sparse_mex, '-^', 'LineWidth', 2);

xlabel('Matrix size n');
ylabel('Time (seconds)');
legend('Gaussian', 'Sparse (MATLAB)', 'Sparse (MEX)');
title('Sketching Performance Scaling');
grid on;


%% testing gaussian vs sparsestack in ID


ns = round(logspace(3, 4.2, 5));  % adjust as needed
k_frac = 0.05;
p = 10;
sketchers = {'gaussian','sparsestack_slow','sparsestack'};

T = zeros(length(ns), length(sketchers));

ntrials = 3;

for ii = 1:length(ns)

    n = ns(ii);
    m = round(n*2.5);
    k = max(10, round(k_frac * n));

    U = orth(randn(n, k));
    V = orth(randn(m, k));
    
    s = exp(-sqrt(0:k-1));   % decaying spectrum
    A = U * diag(s) * V';


    fprintf('n = %d\n', n);

    Anorm = norm(A,'fro');

    for s = 1:length(sketchers)

        t = zeros(ntrials,1);
        err = zeros(ntrials,1);

        for r = 1:ntrials

            tic;
            [Z, rows] = hssutil.inter_decompv3(A, oversample = p, sketcher = sketchers{s});

            t(r) = toc;

            % -----------------------------
            % reconstruction error
            % -----------------------------
            Ahat = Z * A(rows,:);
            err(r) = norm(A - Ahat,'fro') / Anorm;

        end

        T(ii,s) = mean(t);
        E(ii,s) = mean(err);

    end
end

% ---------------- plot ----------------
figure;
loglog(ns, T(:,1), '-o', 'LineWidth', 2); hold on;
loglog(ns, T(:,2), '-s', 'LineWidth', 2);
loglog(ns, T(:,3), '-^', 'LineWidth', 2);

xlabel('n');
ylabel('ID time (s)');
legend(sketchers, 'Location', 'northwest');
title('ID Performance: Gaussian vs SparseStack');
grid on;

figure;
loglog(ns, E(:,1), '-o', 'LineWidth', 2); hold on;
loglog(ns, E(:,2), '-s', 'LineWidth', 2);
loglog(ns, E(:,3), '-^', 'LineWidth', 2);

xlabel('n');
ylabel('relative Frobenius error');
legend(sketchers, 'Location', 'southwest');
title('ID Accuracy: Gaussian vs SparseStack');
grid on;



%% Gaussian vs SparseStack timing 
rng(42);
ns      = round(logspace(3, 4.2, 10));
k_frac  = 0.05;
p       = 10;
ntrials = 7;

sketchers  = {'gaussian', 'sparsestack'};
ns_labels  = {'Gaussian', 'SparseStack'};
T = zeros(length(ns), length(sketchers));
% E = zeros(length(ns), length(sketchers));

for ii = 1:length(ns)
    n    = ns(ii);
    % k    = max(10, round(k_frac * n));
    % U    = orth(randn(n, k));
    % V    = orth(randn(n, k));
    % sv   = exp(-0.2*(0:k-1));          % renamed from s to avoid shadowing
    % A    = U * diag(sv) * V';
    % Cauchy matrix with interlaced points on [-10,10]
    % x = -10 + 20*(0:n-1)'/n;
    % y = -10 + 10/n + 20*(0:n-1)'/n;
    % A = 1 ./ (x - y');           % n x n, full rank, numerically low-rank off-diagonals
    A = randn(n,n);
    Anorm = norm(A, 'fro');
    fprintf('n = %d, k = %d\n', n, k);

    for si = 1:length(sketchers)
        t_vec   = zeros(ntrials, 1);
        err_vec = zeros(ntrials, 1);
        for r = 1:ntrials
            %gaussian
            if si == 1
                tic;
                Omega = randn(n,p) / sqrt(p);
                Y = A * Omega;
                t_vec(r) = toc;
            %sparestack
            elseif si == 2
                zeta = min(8, min(size(A)));
                tic
                S = sparsesign_slow(n,p,zeta);
                Y = A * S;
                t_vec(r) = toc;
            end
            tic;
            % [Z, rows] = hssutil.inter_decompv3(A, 'oversample', p, ...
            %                                'sketcher', sketchers{si});
            % t_vec(r) = toc;
            % Ahat        = Z * A(rows, :);
            % err_vec(r)  = norm(A - Ahat, 'fro') / Anorm;
        end
        % drop best and worst timing trial to reduce noise
        t_vec_sorted = sort(t_vec);
        T(ii, si) = mean(t_vec_sorted(2:end-1));
        % E(ii, si) = mean(err_vec);
    end
end

% ---- publishable timing plot ----
colors = [0.2157 0.4941 0.7216;   % blue  (Gaussian)
          0.8941 0.1020 0.1098];  % red   (SparseStack)
markers = {'o', '^'};

fig1 = figure('Units','inches','Position',[1 1 5.5 3.8]);
ax   = axes('Parent', fig1);
hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
ax.GridAlpha       = 0.2;
ax.GridLineStyle   = '--';
ax.XScale          = 'log';
ax.YScale          = 'log';

for si = 1:length(sketchers)
    loglog(ax, ns, T(:,si), ...
        ['-' markers{si}], ...
        'Color',           colors(si,:), ...
        'LineWidth',       1.8, ...
        'MarkerSize',      6, ...
        'MarkerFaceColor', colors(si,:));
end

% reference line O(n^2)
ref_x = [ns(1) ns(end)];
ref_y = T(1,1) * (ref_x / ns(1)).^2;
loglog(ax, ref_x, ref_y, 'k--', 'LineWidth', 1.2);

xlabel(ax, '$n$',           'Interpreter', 'latex', 'FontSize', 14);
ylabel(ax, 'Time (s)',      'Interpreter', 'latex', 'FontSize', 14);
legend(ax, [ns_labels, {'$O(n^2)$ reference'}], ...
    'Interpreter', 'latex', ...
    'FontSize',    12, ...
    'Location',    'northwest');
ax.FontSize = 12;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XMinorGrid', 'off', 'YMinorGrid', 'off');

%exportgraphics(fig1, 'id_timing_comparison.pdf', 'ContentType', 'vector');

% % ---- accuracy plot (same standards) ----
% fig2 = figure('Units','inches','Position',[1 1 5.5 3.8]);
% ax2  = axes('Parent', fig2);
% hold(ax2, 'on'); box(ax2, 'on'); grid(ax2, 'on');
% ax2.GridAlpha      = 0.2;
% ax2.GridLineStyle  = '--';
% ax2.XScale         = 'log';
% ax2.YScale         = 'log';
% 
% for si = 1:length(sketchers)
%     loglog(ax2, ns, E(:,si), ...
%         ['-' markers{si}], ...
%         'Color',           colors(si,:), ...
%         'LineWidth',       1.8, ...
%         'MarkerSize',      6, ...
%         'MarkerFaceColor', colors(si,:));
% end
% 
% xlabel(ax2, '$n$',                              'Interpreter', 'latex', 'FontSize', 14);
% ylabel(ax2, 'Relative Frobenius error',         'Interpreter', 'latex', 'FontSize', 14);
% legend(ax2, ns_labels, ...
%     'Interpreter', 'latex', ...
%     'FontSize',    12, ...
%     'Location',    'southwest');
% ax2.FontSize = 12;
% ax2.TickLabelInterpreter = 'latex';
% set(ax2, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
% 
% exportgraphics(fig2, 'id_accuracy_comparison.pdf', 'ContentType', 'vector');


%% Gaussian vs SparseStack ID timing and accuracy
rng(42);
ns      = round(logspace(3, 4.2, 10));
k_frac  = 0.05;
p       = 10;
ntrials = 7;

sketchers  = {'gaussian', 'sparsestack'};
ns_labels  = {'Gaussian', 'SparseStack'};
T = zeros(length(ns), length(sketchers));
E = zeros(length(ns), length(sketchers));

for ii = 1:length(ns)
    n    = ns(ii);
    k    = max(10, round(k_frac * n));
    U    = orth(randn(n, k));
    V    = orth(randn(n, k));
    sv   = exp(-0.2*(0:k-1));          
    A    = U * diag(sv) * V';
    Anorm = norm(A, 'fro');
    fprintf('n = %d, k = %d\n', n, k);

    for si = 1:length(sketchers)
        t_vec   = zeros(ntrials, 1);
        err_vec = zeros(ntrials, 1);
        for r = 1:ntrials
            tic;
            [Z, rows] = hssutil.inter_decompv3(A, 'oversample', p, ...
                                           'sketcher', sketchers{si});
            t_vec(r) = toc;
            Ahat        = Z * A(rows, :);
            err_vec(r)  = norm(A - Ahat, 'fro') / Anorm;
        end
        % drop best and worst timing trial to reduce noise
        t_vec_sorted = sort(t_vec);
        T(ii, si) = mean(t_vec_sorted(2:end-1));
        E(ii, si) = mean(err_vec);
    end
end

% ---- publishable timing plot ----
colors = [0.2157 0.4941 0.7216;   % blue  (Gaussian)
          0.8941 0.1020 0.1098];  % red   (SparseStack)
markers = {'o', '^'};

fig1 = figure('Units','inches','Position',[1 1 5.5 3.8]);
ax   = axes('Parent', fig1);
hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
ax.GridAlpha       = 0.2;
ax.GridLineStyle   = '--';
ax.XScale          = 'log';
ax.YScale          = 'log';

for si = 1:length(sketchers)
    loglog(ax, ns, T(:,si), ...
        ['-' markers{si}], ...
        'Color',           colors(si,:), ...
        'LineWidth',       1.8, ...
        'MarkerSize',      6, ...
        'MarkerFaceColor', colors(si,:));
end

% reference line O(n^2)
ref_x = [ns(1) ns(end)];
ref_y = T(1,1) * (ref_x / ns(1)).^2;
loglog(ax, ref_x, ref_y, 'k--', 'LineWidth', 1.2);

xlabel(ax, '$n$',           'Interpreter', 'latex', 'FontSize', 14);
ylabel(ax, 'Time (s)',      'Interpreter', 'latex', 'FontSize', 14);
legend(ax, [ns_labels, {'$O(n^2)$ reference'}], ...
    'Interpreter', 'latex', ...
    'FontSize',    12, ...
    'Location',    'northwest');
ax.FontSize = 12;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XMinorGrid', 'off', 'YMinorGrid', 'off');

exportgraphics(fig1, 'id_timing_comparison.pdf', 'ContentType', 'vector');

% ---- accuracy plot (same standards) ----
fig2 = figure('Units','inches','Position',[1 1 5.5 3.8]);
ax2  = axes('Parent', fig2);
hold(ax2, 'on'); box(ax2, 'on'); grid(ax2, 'on');
ax2.GridAlpha      = 0.2;
ax2.GridLineStyle  = '--';
ax2.XScale         = 'log';
ax2.YScale         = 'log';

for si = 1:length(sketchers)
    loglog(ax2, ns, E(:,si), ...
        ['-' markers{si}], ...
        'Color',           colors(si,:), ...
        'LineWidth',       1.8, ...
        'MarkerSize',      6, ...
        'MarkerFaceColor', colors(si,:));
end

xlabel(ax2, '$n$',                              'Interpreter', 'latex', 'FontSize', 14);
ylabel(ax2, 'Relative Frobenius error',         'Interpreter', 'latex', 'FontSize', 14);
legend(ax2, ns_labels, ...
    'Interpreter', 'latex', ...
    'FontSize',    12, ...
    'Location',    'southwest');
ax2.FontSize = 12;
ax2.TickLabelInterpreter = 'latex';
set(ax2, 'XMinorGrid', 'off', 'YMinorGrid', 'off');

exportgraphics(fig2, 'id_accuracy_comparison.pdf', 'ContentType', 'vector');

%% ID (Gaussian) vs ID (SparseStack) vs SVD — accuracy vs rank k
rng(42);

n       = 3000;                          % fixed matrix size
m       = 7500;
ks      = round(logspace(1, 2.3, 12));   % ranks to sweep  (~10 … 200)
p       = 10;                            % oversampling for ID
ntrials = 7;

methods     = {'gaussian', 'sparsestack', 'svd'};
method_labels = {'ID (Gaussian)', 'ID (SparseStack)', 'SVD'};

% ground-truth matrix: slow singular-value decay so errors are meaningful
% k_true = max(ks);
% U  = orth(randn(n, n));
% V  = orth(randn(m, n));
% %sv = exp(-sqrt(0:k_true-1));
% sv = exp(-(0:n-1).^0.6);
% %sv = 0.99.^(0:n-1);
% A  = U * diag(sv) * V';
% Anorm = norm(A, 'fro');
% x = (1:n)';
% y = (1:m)'+0.5;
% A = 1 ./ (x + y');
[~,H] = rectsampler(n,m,200,400);

A = blockbuilder(H.A11.A12,H.A11);

E = zeros(length(ks), length(methods));
T = zeros(length(ks), length(methods));
for ii = 1:length(ks)
    k = ks(ii);
    fprintf('k = %d\n', k);

    err_vec = zeros(ntrials, 1);
    time_vec = zeros(ntrials, 1);

    % --- ID (Gaussian) ---
    for r = 1:ntrials
        tic;
        [Z, rows] = hssutil.inter_decompv3(A, 'oversample', p, 'sketcher', 'gaussian', 'ctype', 'k', 'cval', k);
        time_vec(r) = toc;
        Ahat = Z * A(rows, :);
        err_vec(r) = norm(A - Ahat, 'fro') / Anorm;
    end
    E(ii, 1) = mean(err_vec);
    T(ii, 1) = mean(time_vec);
    % --- ID (SparseStack) ---
    for r = 1:ntrials
        tic;
        [Z, rows] = hssutil.inter_decompv3(A, 'oversample', p, 'sketcher', 'sparsestack', 'ctype', 'k', 'cval', k);
        time_vec(r) = toc;
        Ahat = Z * A(rows, :);
        err_vec(r) = norm(A - Ahat, 'fro') / Anorm;

    end
    E(ii, 2) = mean(err_vec);
    T(ii, 2) = mean(time_vec);
    % --- Truncated SVD (best possible rank-k approximation) ---
    tic;
    [Us, Ss, Vs] = svds(A, k);
    T(ii, 3) = toc;
    Ahat_svd = Us * Ss * Vs';
    E(ii, 3) = norm(A - Ahat_svd, 'fro') / Anorm;

end

% ---- plot ----
colors = [0.2157 0.4941 0.7216;   % blue  — Gaussian
          0.8941 0.1020 0.1098;   % red   — SparseStack
          0.1725 0.6275 0.1725];  % green — SVD
markers  = {'o', '^', 's'};
lstyles  = {'-', '-', '--'};      % dashed for SVD (it's a bound, not a method)

fig = figure('Units', 'inches', 'Position', [1 1 5.5 3.8]);
ax  = axes('Parent', fig);
hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
ax.GridAlpha     = 0.2;
ax.GridLineStyle = '--';
ax.XScale        = 'log';
ax.YScale        = 'log';

for mi = 1:length(methods)
    loglog(ax, ks, E(:, mi), ...
        [lstyles{mi} markers{mi}], ...
        'Color',           colors(mi, :), ...
        'LineWidth',       1.8, ...
        'MarkerSize',      6, ...
        'MarkerFaceColor', colors(mi, :));
end

xlabel(ax, '$k$  (target rank)',          'Interpreter', 'latex', 'FontSize', 14);
ylabel(ax, 'Relative Frobenius error',    'Interpreter', 'latex', 'FontSize', 14);
legend(ax, method_labels, ...
    'Interpreter', 'latex', ...
    'FontSize',    12, ...
    'Location',    'southwest');
ax.FontSize = 12;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XMinorGrid', 'off', 'YMinorGrid', 'off');

%exportgraphics(fig, 'id_accuracy_vs_rank.pdf', 'ContentType', 'vector');

% ---- timing plot ----
fig2 = figure('Units', 'inches', 'Position', [1 1 5.5 3.8]);
ax2  = axes('Parent', fig2);
hold(ax2, 'on'); box(ax2, 'on'); grid(ax2, 'on');
ax2.GridAlpha = 0.2; ax2.GridLineStyle = '--';
ax2.XScale = 'log'; ax2.YScale = 'log';

for mi = 1:length(methods)
    loglog(ax2, ks, T(:, mi), [lstyles{mi} markers{mi}], ...
        'Color', colors(mi,:), 'LineWidth', 1.8, ...
        'MarkerSize', 6, 'MarkerFaceColor', colors(mi,:));
end
xlabel(ax2, '$k$  (target rank)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel(ax2, 'Time (s)',           'Interpreter', 'latex', 'FontSize', 14);
legend(ax2, method_labels, 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northwest');
ax2.FontSize = 12; ax2.TickLabelInterpreter = 'latex';
set(ax2, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
%exportgraphics(fig2, 'id_timing_vs_rank.pdf', 'ContentType', 'vector');

%% trying to manually create unitary matrices to lower triangularize

M = 30;
N = 60;
k = 5;
blocksize = 10;

[Zfunc, HZ, dH] = rectsampler(M,N,k,blocksize,kernel_type = 'oscillatory');
Z = Zfunc(1:M,1:N);


% [Q11,~] = qr(Z(1:5,1:15));




%% Helper functions


function [Z_handle, HZ, varargout] = rectsampler(M, N, k, blocksize, varargin)
    % Create structured low-rank rectangular matrix
    % Optional: kernel_type = 'random' (default), 'oscillatory', 'decay'
    
    p = inputParser;
    addParameter(p, 'kernel_type', 'random', @ischar);
    addParameter(p, 'freq', 10, @isnumeric);  % for oscillatory
    addParameter(p, 'decay_rate', 0.1, @isnumeric);  % for decay
    parse(p, varargin{:});
    
    kernel_type = p.Results.kernel_type;
    
    % Create geometry (like your Hankel example)
    t_row = linspace(0, 1, M+1); t_row(end) = [];
    t_col = linspace(0, 1, N+1); t_col(end) = [];
    
    % Store geometry for block function
    geom.t_row = t_row(:);
    geom.t_col = t_col(:);
    geom.k = k;
    geom.kernel_type = kernel_type;
    geom.freq = p.Results.freq;
    geom.decay_rate = p.Results.decay_rate;
    
    % Create function handle
    Z_handle = @(I, J) rect_kernel_block(I, J, geom);
    
    % Build HSS from function handle
    tic;
    HZ = hss(Z_handle, blocksize=blocksize, sizeA=[M, N]);
    timing = toc;
    
    % Apply diagonal modification
    [HZ,dH] = diagmodify(HZ);
    
    Z_handle = @(I,J) rect_kernel_block(I,J,geom) + dH(I,J);

    % Return timing if requested
    if nargout == 3
        varargout{1} = dH;
    elseif nargout == 4
        varargout{1} = dH;
        varargout{2} = timing;
    end
end

function B = rect_kernel_block(I, J, geom)
    % Compute block of low-rank structured matrix
    
    ti = geom.t_row(I);
    tj = geom.t_col(J);
    k = geom.k;
    
    switch geom.kernel_type
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