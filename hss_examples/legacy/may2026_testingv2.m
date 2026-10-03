%% =========================================================================
%%  hss_ulv_recursive_test.m
%%  ─────────────────────────────────────────────────────────────────────────
%%  Applies the HSS-ULV unitary transforms ONE NODE AT A TIME through the
%%  tree.  Two hard rules:
%%
%%   (1) Each Q touches ONLY its own node's column block.
%%       The full N×N working matrix Zw is never passed to qr().
%%
%%   (2) Permutations are NEVER FORMED as matrices.
%%       Visualisation reorderings are plain index vectors used as
%%       column subscripts:  ZvisN = Zw(:, permN)
%%
%%  After Phase 1 (leaf LQ) + Phase 2 (level-1 LQ) the result has the
%%  hierarchical ULV structure needed for the recursive min-norm solve:
%%
%%    active block (M×M) : block lower-triangular  +  rank-k cross-off-diag
%%    null   block (M×(N-M)) : low-rank (rank ≤ k per row-block)
%%
%%  No dense root QR (Phase 3) is applied.  The cross-off-diagonal is the
%%  Schur complement connection resolved by the recursive solve.
%% =========================================================================

clear; close all; clc

%% ── 0.  Build test problem ───────────────────────────────────────────────
M = 30;  N = 60;  k = 5;  blocksize = 10;
[Zfunc, HZ, ~] = rectsampler(M, N, k, blocksize, 'kernel_type', 'oscillatory');
Z  = Zfunc(1:M, 1:N);
Zw = Z;          % working copy  – only column-block transforms applied
V  = eye(N);     % accumulated right unitary  (optional; for verification)

leaves = gather_leaves(HZ);
nL     = numel(leaves);
fprintf('Leaves: %d\n', nL);

%% ═════════════════════════════════════════════════════════════════════════
%%  PHASE 1 – LEAF LQ
%%
%%  For leaf i with diagonal block D_i = Zw(ri, ci)  (m_i × n_i, m_i < n_i):
%%
%%    [Q_i, ~] = qr( D_i' )          ← QR of the m_i×n_i block transposed
%%    Zw(:, ci) = Zw(:, ci) * Q_i    ← column-block transform; only ci touched
%%
%%  Result inside the block:   Zw(ri, ci(1:m_i))     = L_i  (lower tri)
%%                             Zw(ri, ci(m_i+1:end)) ≈ 0    (null-space of D_i)
%%  Outside (off-diagonal rows): non-zero, rank ≤ k  (generator structure)
%% ═════════════════════════════════════════════════════════════════════════
fprintf('\n── Phase 1: Leaf LQ ─────────────────────────────────────────\n');

act = [];   % active col indices  (span row-space of each D_i)
nul = [];   % null   col indices  (span null-space of each D_i)

for i = 1:nL
    ri = leaves(i).rows;   mi = numel(ri);
    ci = leaves(i).cols;   ni = numel(ci);

    % QR of D_i'  (n_i × m_i) – only the m_i × n_i leaf block is read
    [Qi, ~]   = qr( Zw(ri, ci)' );

    % Right transform applied to column block ci ONLY
    Zw(:, ci) = Zw(:, ci) * Qi;
    V(:,  ci) = V(:,  ci) * Qi;

    act = [act, ci(1:mi)    ];
    nul = [nul, ci(mi+1:ni) ];

    fprintf('  Leaf %d  rows[%d:%d]  cols[%d:%d]  D:%d×%d  active:%d  null:%d\n', ...
            i, ri(1),ri(end), ci(1),ci(end), mi,ni, mi, ni-mi);
end

% ── Snapshot 1: column index vector  [active (M cols) | null (N-M cols)] ──
% perm1 is a plain index vector, never a matrix.
perm1 = [act, nul];
Zvis1 = Zw(:, perm1);   % Zw reindexed by subscript – no permutation matrix

% Verify: diagonal null columns should be zero
diag_null_res = 0;
for i = 1:nL
    ri = leaves(i).rows; ci = leaves(i).cols; mi = numel(ri);
    diag_null_res = diag_null_res + norm(Zw(ri, ci(mi+1:end)), 'fro');
end
fprintf('\n  ||diag blocks at null cols||_F = %.2e  (expect ~0)\n', diag_null_res);


%% ═════════════════════════════════════════════════════════════════════════
%%  PHASE 2 – LEVEL-1 LQ
%%
%%  For level-1 node j (HZ.A11 or HZ.A22):
%%    ca  = active cols from Phase 1 that fall inside node j's column range
%%    The sub-block Zw(ri_j, ca) is SQUARE  (|ri_j| × |ri_j| = M/2 × M/2)
%%    because the leaf LQ reduced each leaf's n_i cols to m_i active cols,
%%    and  Σ m_i = M/2  within each level-1 subtree.
%%
%%    [Q_j, ~] = qr( Zw(ri_j, ca)' )   ← QR of the M/2 × M/2 block
%%    Zw(:, ca) = Zw(:, ca) * Q_j       ← only columns ca are touched
%% ═════════════════════════════════════════════════════════════════════════
fprintf('\n── Phase 2: Level-1 LQ ──────────────────────────────────────\n');
nd_all   = {HZ.A11, HZ.A22};
for j = 1:2
    nd   = nd_all{j};

    ri_j = nd.Ir(1):nd.Ir(2);
    ci_j = nd.Ic(1):nd.Ic(2);

    ca = act( ismember(act, ci_j) );   % active cols inside this node

    % QR of the square sub-block Zw(ri_j, ca)  (M/2 × M/2)
    [Qj, ~]   = qr( Zw(ri_j, ca)' );

    % Right transform applied to column block ca ONLY
    Zw(:, ca) = Zw(:, ca) * Qj;
    V(:,  ca) = V(:,  ca) * Qj;

    Bj = Zw(ri_j, ca);
    ut = norm(triu(Bj,1),'fro') / (norm(Bj,'fro') + eps);
    fprintf('  Node %d  rows[%d:%d]  |ca|=%d  B:%d×%d  ||upper-tri||/||B||=%.2e\n', ...
            j, ri_j(1),ri_j(end), numel(ca), numel(ri_j),numel(ca), ut);
end

% ── Snapshot 2: [A11-active | A22-active | null] ──────────────────────────
% Again: a plain index vector, never a matrix.
act11 = act( ismember(act, HZ.A11.Ic(1):HZ.A11.Ic(2)) );
act22 = act( ismember(act, HZ.A22.Ic(1):HZ.A22.Ic(2)) );
M2    = numel(act11);   % = M/2 = 15

perm2 = [act11, act22, nul];
Zvis2 = Zw(:, perm2);


%% ═════════════════════════════════════════════════════════════════════════
%%  DIAGNOSTICS – structure after Phase 1 + Phase 2 (no Phase 3)
%% ═════════════════════════════════════════════════════════════════════════
fprintf('\n── Diagnostics (no Phase 3 applied) ────────────────────────\n');

B11 = Zw(1:M2,   act11);   % M/2 × M/2 – should be lower triangular
B22 = Zw(M2+1:M, act22);   % M/2 × M/2 – should be lower triangular
G12 = Zw(1:M2,   act22);   % M/2 × M/2 – cross-off-diag, rank ≤ k
G21 = Zw(M2+1:M, act11);   % M/2 × M/2 – cross-off-diag, rank ≤ k

fprintf('  Diagonal blocks (should be lower triangular):\n');
fprintf('    ||upper-tri(B11)||/||B11|| = %.2e\n', ...
        norm(triu(B11,1),'fro')/(norm(B11,'fro')+eps));
fprintf('    ||upper-tri(B22)||/||B22|| = %.2e\n', ...
        norm(triu(B22,1),'fro')/(norm(B22,'fro')+eps));

sv12 = svd(G12);
sv21 = svd(G21);
fprintf('\n  Cross-off-diagonal singular values (rank ≤ k=%d expected):\n', k);
fprintf('    G12:  '); fprintf('%.2e ', sv12(1:min(k+2,end))'); fprintf('\n');
fprintf('    G21:  '); fprintf('%.2e ', sv21(1:min(k+2,end))'); fprintf('\n');
fprintf('    G12 effective rank at tol 1e-10: %d\n', sum(sv12 > 1e-10));
fprintf('    G21 effective rank at tol 1e-10: %d\n', sum(sv21 > 1e-10));


%% ═════════════════════════════════════════════════════════════════════════
%%  VISUALISATION
%%  Six panels:
%%   [1] Original Z
%%   [2] After Phase 1  (index-permuted view: act | null)
%%   [3] After Phase 2  (index-permuted view: act11 | act22 | null)
%%   [4] Active block M×M  (block lower-tri + rank-k cross)
%%   [5] Null block  M×(N-M)  (low-rank per row-block)
%%   [6] Singular value decay of cross-off-diagonal G12 and G21
%% ═════════════════════════════════════════════════════════════════════════
figure('Name','Recursive HSS-ULV (Phase 1 + Phase 2, no dense Phase 3)', ...
       'Units','normalized','Position',[.01 .05 .96 .87]);

%% helper: draw a vertical white dashed line (compatible, no xline needed)
vline = @(ax,x) plot(ax, [x x], [.5 M+.5], 'w--', 'LineWidth', 1.5);
hline = @(ax,y) plot(ax, [.5 N+.5], [y y], 'w--', 'LineWidth', 1.5);

%% Panel 1 – original
ax1 = subplot(2,3,1);
imagesc(abs(Z)); colorbar; colormap(ax1,'parula');
title('Step 0 – Z  (original)','FontWeight','bold','FontSize',9);
xlabel('col'); ylabel('row');

%% Panel 2 – after Phase 1  (columns: active | null)
ax2 = subplot(2,3,2);
imagesc(abs(Zvis1)); colorbar; colormap(ax2,'parula'); hold on;
vline(ax2, M+.5);
title({'Step 1 – after Leaf LQ', 'cols: [active (1:M) | null (M+1:N)]'}, ...
      'FontWeight','bold','FontSize',9);
xlabel('col (permuted)'); ylabel('row');

%% Panel 3 – after Phase 2  (columns: A11-active | A22-active | null)
ax3 = subplot(2,3,3);
imagesc(abs(Zvis2)); colorbar; colormap(ax3,'parula'); hold on;
vline(ax3, M2+.5);
vline(ax3, M+.5);
hline(ax3, M2+.5);
title({'Step 2 – after Level-1 LQ', 'cols: [A11-act | A22-act | null]'}, ...
      'FontWeight','bold','FontSize',9);
xlabel('col (permuted)'); ylabel('row');

%% Panel 4 – active block zoom  (the M×M square that the solve works with)
ax4 = subplot(2,3,4);
imagesc(abs(Zvis2(:, 1:M))); colorbar; colormap(ax4,'parula'); hold on;
vline(ax4, M2+.5);
hline(ax4, M2+.5);
title({'Active block  Z_w(:, act)  M×M', ...
       '[L_{11}  G_{12}; G_{21}  L_{22}]   rank(G_{12,21}) ≤ k'}, ...
      'FontWeight','bold','FontSize',9);
xlabel('active col'); ylabel('row');

%% Panel 5 – null block  (M × (N-M))
ax5 = subplot(2,3,5);
imagesc(abs(Zvis2(:, M+1:end))); colorbar; colormap(ax5,'parula'); hold on;
hline(ax5, M2+.5);
title({'Null block  Z_w(:, null)  M×(N-M)', 'rank ≤ k per row-block'}, ...
      'FontWeight','bold','FontSize',9);
xlabel('null col'); ylabel('row');

%% Panel 6 – singular value decay of cross-off-diagonal
ax6 = subplot(2,3,6);
semilogy(sv12, 'bo-', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','G_{12}'); hold on;
semilogy(sv21, 'rs-', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','G_{21}');
plot([k+.5 k+.5], ylim, 'k--', 'LineWidth',1);
text(k+.7, max(sv12)*0.5, sprintf('k=%d',k), 'FontSize',9);
legend('Location','southwest'); grid on;
title({'Singular values of G_{12}, G_{21}', ...
       '(cross-off-diagonal; this is the Schur complement)'}, ...
      'FontWeight','bold','FontSize',9);
xlabel('index'); ylabel('\sigma_i');

sgtitle(sprintf(...
  'HSS-ULV  [M=%d  N=%d  k=%d]  –  recursive tree transforms, no dense Phase 3', ...
  M, N, k), 'FontSize',11,'FontWeight','bold');


%% ═════════════════════════════════════════════════════════════════════════
%%  LOCAL FUNCTIONS
%% ═════════════════════════════════════════════════════════════════════════

function lv = gather_leaves(node)
%GATHER_LEAVES  Collect HSS leaf nodes left-to-right.
    lv = struct('rows',{},'cols',{});
    lv = rec(node, lv);
end

function lv = rec(nd, lv)
    if nd.isleaf
        lv(end+1).rows = nd.Ir(1):nd.Ir(2);
        lv(end  ).cols = nd.Ic(1):nd.Ic(2);
    else
        lv = rec(nd.A11, lv);
        lv = rec(nd.A22, lv);
    end
end


%% helper functions
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