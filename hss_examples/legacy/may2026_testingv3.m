%% =========================================================================
%%  hss_ulv_test.m
%%  ─────────────────────────────────────────────────────────────────────────
%%  At each leaf pair within a level-1 subtree:
%%
%%   LEFT   QL of the off-diagonal ROW GENERATOR  nd.A12.Z  (m×k)
%%          → zeros the top (m−k) rows of that generator
%%          → by HSS nesting, those same rows see zero in EVERY off-diagonal
%%            so they are fully decoupled
%%
%%   RIGHT  full LQ of the (now-Q'-transformed) dense diagonal block
%%          → diagonal block becomes lower triangular
%%          → the decoupled rows can be solved immediately (bsolve)
%%
%%   RECOMPRESS  the remaining k rows per leaf are the Schur complement;
%%               collect them into a 4k×4k system and solve in one shot
%%
%%  Nothing full-matrix is ever formed.  Permutations are index vectors.
%% =========================================================================
clear; close all; clc

%% ── 0.  Setup ────────────────────────────────────────────────────────────
M = 30;  N = 60;  k = 5;  blocksize = 10;
[Zfunc, HZ, ~] = rectsampler(M, N, k, blocksize, 'kernel_type', 'oscillatory');
Z  = Zfunc(1:M, 1:N);
rng(0);  b = randn(M, 1);

Zw     = Z;
b_work = b;
wbar   = zeros(M, 1);
wind   = false(M, 1);
Ss      = {};   % Ss{i}      = S_i (n_i×m_i)  right factor for leaf i
CRanges = {};   % CRanges{i} = col indices of leaf i in full system

%% ── Leaf bookkeeping ─────────────────────────────────────────────────────
leaves = gather_leaves(HZ);
nL     = numel(leaves);
m_vec  = arrayfun(@(i) numel(leaves(i).rows), 1:nL);
n_vec  = arrayfun(@(i) numel(leaves(i).cols), 1:nL);
ws     = [0, cumsum(m_vec(1:end-1))];        % wbar start indices (0-based)

ri = {leaves.rows};
ci = {leaves.cols};

fprintf('M=%d  N=%d  k=%d   leaves=%d   Schur dim=%d\n', M,N,k,nL,nL*k);

%% ═════════════════════════════════════════════════════════════════════════
%%  LEAF-PAIR LOOP
%%  One pass per level-1 subtree.  The loop body is the same each time;
%%  only the generator source (HZ.A11 vs HZ.A22) and leaf indices change.
%% ═════════════════════════════════════════════════════════════════════════
level1_nodes = {HZ.A11, HZ.A22};
leaf_pairs   = [1 2; 3 4];            % which leaves live in each subtree

for p = 1:2
    nd = level1_nodes{p};
    ia = leaf_pairs(p,1);   ib = leaf_pairs(p,2);

    fprintf('\n── Subtree %d  (leaves %d, %d) ─────────────────────────────\n',p,ia,ib);

    % ── LEFT: QL of the off-diagonal row generators ──────────────────────
    %
    %  nd.A12.Z  is the m_ia×k row generator:
    %    "how leaf ia's rows contribute to the off-diagonal connecting to leaf ib"
    %  QL of this k-column matrix zeros out its top (m_ia − k) rows.
    %  Because every off-diagonal block that touches leaf ia's rows inherits
    %  this generator (HSS nested basis), those (m_ia − k) rows become
    %  decoupled from ALL off-diagonal blocks simultaneously.
    %
    %  nd.A21.Z  plays the same role for leaf ib.

    Qa = leafcompressQ(nd.A12.Z);          % m_ia × m_ia  (QL of m_ia×k matrix)
    Zw(ri{ia}, :)  = Qa' * Zw(ri{ia}, :); % row-block transform (touches ri{ia} only)
    b_work(ri{ia}) = Qa' * b_work(ri{ia});

    Qb = leafcompressQ(nd.A21.Z);          % m_ib × m_ib
    Zw(ri{ib}, :)  = Qb' * Zw(ri{ib}, :);
    b_work(ri{ib}) = Qb' * b_work(ri{ib});

    % Verify decoupling: off-diagonal block at the top rows should be ≈ 0
    t_ia = m_vec(ia) - k;
    t_ib = m_vec(ib) - k;
    fprintf('  Left QL from generators:\n');
    fprintf('    leaf%d: ||Zw(ri_a(1:%d), ci_b)||_F = %.2e  (expect ~0)\n', ...
            ia, t_ia, norm(Zw(ri{ia}(1:t_ia), ci{ib})));
    fprintf('    leaf%d: ||Zw(ri_b(1:%d), ci_a)||_F = %.2e  (expect ~0)\n', ...
            ib, t_ib, norm(Zw(ri{ib}(1:t_ib), ci{ia})));

    % ── RIGHT: full LQ of each (Q'-transformed) diagonal block ────────────
    %
    %  Zw(ri{ia}, ci{ia}) is now  Qa'·D_ia  (m_ia × n_ia, m_ia < n_ia).
    %  Full QR of its transpose gives D_ia_new · Q_right = [L_ia | 0]
    %  where L_ia is m_ia×m_ia lower triangular.
    %  Store S_ia = Q_right(:,1:m_ia) (n_ia×m_ia) for min-norm assembly.

    [Qa_r, ~]    = qr(Zw(ri{ia}, ci{ia})');   % full QR of (n_ia × m_ia)
    Zw(:, ci{ia}) = Zw(:, ci{ia}) * Qa_r;     % col-block transform (ci{ia} only)
    Ss{end+1}      = Qa_r(:, 1:m_vec(ia));    % n_ia × m_ia
    CRanges{end+1} = ci{ia};

    [Qb_r, ~]    = qr(Zw(ri{ib}, ci{ib})');
    Zw(:, ci{ib}) = Zw(:, ci{ib}) * Qb_r;
    Ss{end+1}      = Qb_r(:, 1:m_vec(ib));
    CRanges{end+1} = ci{ib};

    D_ia = Zw(ri{ia}, ci{ia}(1:m_vec(ia)));   % m_ia × m_ia   (should be lower tri)
    D_ib = Zw(ri{ib}, ci{ib}(1:m_vec(ib)));
    fprintf('  Right LQ on diagonal blocks:\n');
    fprintf('    leaf%d: ||upper-tri(D)||/||D|| = %.2e\n', ia, ...
            norm(triu(D_ia,1),'fro')/(norm(D_ia,'fro')+eps));
    fprintf('    leaf%d: ||upper-tri(D)||/||D|| = %.2e\n', ib, ...
            norm(triu(D_ib,1),'fro')/(norm(D_ib,'fro')+eps));

    % ── BSOLVE: solve the decoupled rows immediately ──────────────────────
    %
    %  After left + right transforms, row ri{ia}(1:t_ia) satisfies:
    %    D_ia(1:t_ia, 1:t_ia) · wbar = b_work(ri{ia}(1:t_ia))   EXACTLY
    %  because those rows see zero in every off-diagonal direction.
    %  Solve by forward substitution on the lower-triangular block.

    if t_ia > 0
        wbar(ws(ia)+(1:t_ia)) = D_ia(1:t_ia,1:t_ia) \ b_work(ri{ia}(1:t_ia));
        wind(ws(ia)+(1:t_ia)) = true;
    end
    if t_ib > 0
        wbar(ws(ib)+(1:t_ib)) = D_ib(1:t_ib,1:t_ib) \ b_work(ri{ib}(1:t_ib));
        wind(ws(ib)+(1:t_ib)) = true;
    end

    fprintf('  Bsolve: solved %d (leaf%d) + %d (leaf%d) entries\n', t_ia,ia,t_ib,ib);
end

fprintf('\nLeaf level complete:  %d solved   %d → Schur\n', sum(wind), sum(~wind));

%% ═════════════════════════════════════════════════════════════════════════
%%  RECOMPRESS  +  MOVE UP ONE LEVEL
%%  Each leaf contributes k "Schur rows" (the bottom k rows of its ri block)
%%  and k "Schur-active columns" (the last k of its m_i active columns).
%%  Together they form the  4k × 4k  square Schur complement system.
%%
%%  Schur RHS: subtract the already-solved wbar's contribution from the
%%  Schur rows.  wbar is zero at Schur positions so only solved entries act.
%% ═════════════════════════════════════════════════════════════════════════
schur_rows     = [];
schur_act_cols = [];
schur_wbar_idx = [];
act_all        = [];

for i = 1:nL
    ti = m_vec(i) - k;
    schur_rows     = [schur_rows;     ri{i}(ti+1:end)'  ];   % last k rows
    schur_act_cols = [schur_act_cols; ci{i}(ti+1:m_vec(i))'];% last k active cols
    schur_wbar_idx = [schur_wbar_idx; ws(i)+(ti+1:m_vec(i))];
    act_all        = [act_all;        ci{i}(1:m_vec(i))' ];
end

b_schur  = b_work(schur_rows) - Zw(schur_rows, act_all) * wbar;
Zw_schur = Zw(schur_rows, schur_act_cols);   % 4k × 4k

fprintf('\nSchur system: %d×%d   rank=%d   cond=%.2e\n', ...
        size(Zw_schur,1), size(Zw_schur,2), rank(Zw_schur), cond(Zw_schur));

wbar(schur_wbar_idx) = Zw_schur \ b_schur;

%% ═════════════════════════════════════════════════════════════════════════
%%  ASSEMBLY:  x(ci_i) = S_i · wbar_i
%%  S_i (n_i × m_i) maps the m_i compressed active coordinates back to
%%  the n_i original column directions.  Null-space entries stay zero.
%% ═════════════════════════════════════════════════════════════════════════
x = zeros(N, 1);
w_ptr = 0;
for i = 1:numel(Ss)
    mi = size(Ss{i}, 2);
    x(CRanges{i}) = Ss{i} * wbar(w_ptr + (1:mi));
    w_ptr = w_ptr + mi;
end

%% ═════════════════════════════════════════════════════════════════════════
%%  VERIFICATION
%% ═════════════════════════════════════════════════════════════════════════
x_ref = pinv(Z) * b;
fprintf('\n══════════════════════════════════════════════\n');
fprintf('Residual   ||Zx - b||   = %.3e\n', norm(Z*x - b));
fprintf('           ||x||        = %.6f\n', norm(x));
fprintf('  pinv ref ||x_ref||    = %.6f\n', norm(x_ref));
fprintf('  error    ||x-x_ref||  = %.3e\n', norm(x - x_ref));
fprintf('══════════════════════════════════════════════\n');

%% ═════════════════════════════════════════════════════════════════════════
%%  VISUALISATION  (index permutations only – no permutation matrices)
%% ═════════════════════════════════════════════════════════════════════════
act_vis = [];  nul_vis = [];
for i = 1:nL
    act_vis = [act_vis; ci{i}(1:m_vec(i))'      ];
    nul_vis = [nul_vis; ci{i}(m_vec(i)+1:end)'   ];
end
perm_act_nul = [act_vis; nul_vis];    % [active (M cols) | null (N-M cols)]
M2 = sum(m_vec(1:2));                 % M/2 = 15, the A11 subtree active dim

%% Snapshots of Zw at key stages (re-run transforms in order; save here)
% (Zw is already in the fully-transformed state after the loop above.)
Zw_final = Zw(:, perm_act_nul);      % permuted for display

%% Singular values of cross-off-diagonal (shows Schur coupling rank)
G12 = Zw(1:M2,   act_vis(M2+1:M));
G21 = Zw(M2+1:M, act_vis(1:M2));
sv12 = svd(G12);
sv21 = svd(G21);

figure('Name','HSS-ULV  –  left QL from generators, right LQ on blocks', ...
       'Units','normalized','Position',[.02 .05 .94 .87]);

wl = @(ax,x) plot(ax,[x x],[.5 M+.5],'w--','LineWidth',1.5);
hl = @(ax,y) plot(ax,[.5 N+.5],[y y],'w--','LineWidth',1.5);
wa = @(ax,x) plot(ax,[x x],[.5 M+.5],'c--','LineWidth',1.2);
ha = @(ax,y) plot(ax,[.5 M+.5],[y y],'c--','LineWidth',1.2);

%% Panel 1 – original Z
ax1 = subplot(2,3,1); imagesc(abs(Z)); colorbar; colormap(ax1,'parula');
title('Z  (original)','FontWeight','bold','FontSize',9); xlabel('col'); ylabel('row');

%% Panel 2 – Zw after all transforms, permuted [active | null]
ax2 = subplot(2,3,2); imagesc(abs(Zw_final)); colorbar; colormap(ax2,'parula'); hold on;
wl(ax2, M+.5);         % active | null boundary
wa(ax2, M2+.5);        % A11 | A22 within active
ha(ax2, M2+.5);
title({'Zw after leaf-level QL+LQ','cols: [active | null]'},...
      'FontWeight','bold','FontSize',9); xlabel('col (permuted)'); ylabel('row');

%% Panel 3 – active block M×M  (block lower-tri + rank-k Schur coupling)
ax3 = subplot(2,3,3); imagesc(abs(Zw_final(:,1:M))); colorbar; colormap(ax3,'parula'); hold on;
wa(ax3, M2+.5); ha(ax3, M2+.5);
title({'Active block  M×M', '[L_{11}  G_{12}; G_{21}  L_{22}]   rank(G)≤k'},...
      'FontWeight','bold','FontSize',9); xlabel('active col'); ylabel('row');

%% Panel 4 – 4k×4k Schur system (the "next level up")
ax4 = subplot(2,3,4); imagesc(abs(Zw_schur)); colorbar; colormap(ax4,'parula');
title(sprintf('Schur system  %d×%d  (rank %d)',4*k,4*k,rank(Zw_schur)),...
      'FontWeight','bold','FontSize',9); xlabel('Schur col'); ylabel('Schur row');

%% Panel 5 – singular value decay of cross-off-diagonal
ax5 = subplot(2,3,5);
semilogy(sv12,'bo-','LineWidth',1.5,'MarkerSize',5,'DisplayName','G_{12}'); hold on;
semilogy(sv21,'rs-','LineWidth',1.5,'MarkerSize',5,'DisplayName','G_{21}');
plot([k+.5 k+.5],[1e-14 max(sv12)*2],'k--','LineWidth',1);
text(k+.7,sv12(1)*0.5,sprintf('k=%d',k),'FontSize',9);
legend('Location','southwest'); grid on;
title({'Cross-off-diagonal singular values','(Schur coupling; rank ≤ k)'},...
      'FontWeight','bold','FontSize',9);
xlabel('index'); ylabel('\sigma_i');

%% Panel 6 – solution
ax6 = subplot(2,3,6);
plot(x,   'b-', 'LineWidth',1.5,'DisplayName','x_{ULV}'); hold on;
plot(x_ref,'r--','LineWidth',1.2,'DisplayName','pinv(Z)*b');
legend; grid on;
title(sprintf('Solution  ||x-x_{ref}||=%.1e',norm(x-x_ref)),...
      'FontWeight','bold','FontSize',9);
xlabel('component'); ylabel('value');

sgtitle(sprintf(['HSS-ULV  M=%d  N=%d  k=%d  –  ' ...
    'left QL from Z generators · right LQ on D · bsolve as we compress'], M,N,k),...
    'FontSize',11,'FontWeight','bold');

%% ═════════════════════════════════════════════════════════════════════════
%%  LOCAL FUNCTIONS
%% ═════════════════════════════════════════════════════════════════════════

function Q = leafcompressQ(A)
%LEAFCOMPRESSQ  QL via flip-QR:  Q'·A has its top (m−rank(A)) rows ≈ 0.
%  A is m×k (k < m).  Q is m×m unitary.
    [m, n] = size(A);
    [Q, ~] = qr(A * flip(eye(n)));
    Q      = Q * flip(eye(m));
end

function lv = gather_leaves(node)
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