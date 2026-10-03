%% =========================================================================
%%  hss_ulv_minnorm.m
%%  -----------------------------------------------------------------------
%%  ULV factorization of a rectangular HSS matrix for the min-norm problem.
%%
%%  Builds  V (N×N unitary)  and  L (M×M lower triangular)  satisfying
%%
%%              Z · V  =  [ L  |  0_{M×(N–M)} ]
%%
%%  Min-norm solution to  Z·x = b:
%%              x_min  =  V(:,1:M) · (L \ b)
%%
%%  (U = I throughout because we perform only right-unitary transforms.)
%%
%%  ─────────────────────────────────────────────────────────────────────
%%  ALGORITHM  — three bottom-up phases
%%  ─────────────────────────────────────────────────────────────────────
%%  Phase 1  LEAF LQ
%%    For each leaf diagonal block D_i (m_i×n_i, m_i < n_i):
%%      [Q_i,~] = qr(D_i')  →  D_i·Q_i = [L_i | 0]
%%    Apply Q_i to full column block c_i.  First m_i new cols are
%%    "active" (row-space of D_i); remaining n_i–m_i are "null-space".
%%
%%  Phase 2  LEVEL-1 LQ
%%    For each level-1 node (HZ.A11, HZ.A22), collect its Phase-1
%%    active cols.  The sub-block is SQUARE (M/2 × M/2) because
%%    Σ m_i = M/2 per subtree.  Apply LQ to make it lower triangular.
%%    No new null-space columns are produced.
%%
%%  Phase 3  ROOT / FINAL LQ
%%    One dense QR on the full working matrix absorbs the residual
%%    rank-k cross-level blocks and guarantees Z·V = [L|0] exactly.
%%
%%  ─────────────────────────────────────────────────────────────────────
%%  VISUALIZATION PERMUTATIONS  (never stored in V)
%%  ─────────────────────────────────────────────────────────────────────
%%  perm1  after Phase 1:  cols → [active (1:M) | null (M+1:N)]
%%  perm2  after Phase 2:  cols → [A11-act | A22-act | null]
%%  These let you watch the matrix shift toward lower-triangular form.
%%  Build the corresponding permutation matrices with:
%%      P1 = eye(N); P1 = P1(:, perm1);   % Z_work·P1 = Zv1
%%      P2 = eye(N); P2 = P2(:, perm2);   % Z_work·P2 = Zv2
%% =========================================================================

clear; close all; clc

%% ── 0. Build problem ─────────────────────────────────────────────────────
M = 35;  N = 65;  k = 5;  blocksize = 10;
[Zfunc, HZ, ~] = rectsampler(M,N,k,blocksize,'kernel_type','oscillatory');
Z = Zfunc(1:M,1:N);
[M, N] = size(Z);

%% ── Initialise factorisation accumulator ─────────────────────────────────
V  = eye(N);      % right-unitary accumulator;  Z·V → [L|0]
Zw = Z;           % working matrix;              updated in-place each phase

%% ── Collect leaf nodes (left-to-right, i.e. A11 before A22) ─────────────
leaves = gather_leaves(HZ);
nL     = numel(leaves);

fprintf('HSS tree: %d leaf nodes\n', nL)
for i = 1:nL
    fprintf('  Leaf %d  rows[%d:%d]  cols[%d:%d]  D: %d×%d\n', i, ...
        leaves(i).rows(1), leaves(i).rows(end), ...
        leaves(i).cols(1), leaves(i).cols(end), ...
        numel(leaves(i).rows), numel(leaves(i).cols))
end

%% =========================================================================
%%  PHASE 1 – LEAF LQ
%% =========================================================================
fprintf('\n─── Phase 1: Leaf LQ ───────────────────────────────────────────\n')

act = [];   % active    column indices (row-space of each leaf D_i)
nul = [];   % null-space column indices

for i = 1:nL
    ri = leaves(i).rows;  mi = numel(ri);   % m_i rows in this leaf
    ci = leaves(i).cols;  ni = numel(ci);   % n_i cols in this leaf

    % ── Full QR of D_i'  (n_i × m_i matrix)
    %    D_i'  =  Q_i · R_i   →   D_i · Q_i  =  R_i'  =  [L_i | 0]
    %    Q_i : n_i×n_i unitary,   R_i : n_i×m_i upper-tri
    [Qi, ~]   = qr( Zw(ri, ci)' );    % full QR (no 0/econ flag)

    % ── Apply right transform to entire column block c_i
    Zw(:, ci) = Zw(:, ci) * Qi;       % transforms every row
    V(:,  ci) = V(:,  ci) * Qi;       % accumulate into V

    % ── Track active / null-space columns
    act = [act, ci(1:mi)    ];         % cols 1..m_i  → row-space of D_i
    nul = [nul, ci(mi+1:ni) ];         % cols m_i+1.. → null-space of D_i

    % Sanity: diagonal block should now be [L_i | ~0]
    d_null = norm(Zw(ri, ci(mi+1:end)), 'fro');
    fprintf('  Leaf %d: D=%d×%d  act=%d  nul=%d  ||D·null||_F=%.2e\n', ...
            i, mi, ni, mi, ni-mi, d_null)
end

%% ── Visualisation permutation 1: [active (M cols) | null (N-M cols)]
%%    Not applied to V.  Build P1 = eye(N); P1 = P1(:,perm1) if desired.
perm1 = [act, nul];
Zv1   = Zw(:, perm1);   % snapshot for plotting

%% =========================================================================
%%  PHASE 2 – LEVEL-1 LQ
%% =========================================================================
fprintf('\n─── Phase 2: Level-1 LQ ────────────────────────────────────────\n')

lvl1_nodes = {HZ.A11, HZ.A22};

for j = 1:2
    nd   = lvl1_nodes{j};
    ri_j = nd.Ir(1):nd.Ir(2);            % row range of this level-1 node
    ci_j = nd.Ic(1):nd.Ic(2);            % col range of this level-1 node

    % Active cols that belong to this node's column range
    ca = act(ismember(act, ci_j));        % |ca| = |ri_j| = M/2 → SQUARE

    % ── LQ of the square sub-block
    %    B_j · Q_j  =  L_j   (square lower triangular)
    Bj          = Zw(ri_j, ca);
    [Qj, ~]     = qr(Bj');               % full QR of B_j' (square input)
    Zw(:, ca)   = Zw(:, ca) * Qj;
    V(:,  ca)   = V(:,  ca) * Qj;

    % Verify: diagonal block should now be lower triangular
    Bjnew = Zw(ri_j, ca);
    ut    = norm(triu(Bjnew,1),'fro') / (norm(Bjnew,'fro') + eps);
    fprintf('  Node %d  rows[%d:%d]  |ca|=%d  ||upper-tri(B)||/||B||=%.2e\n', ...
            j, ri_j(1), ri_j(end), numel(ca), ut)
end

%% ── Visualisation permutation 2: [A11-active | A22-active | null]
act11 = act(ismember(act, HZ.A11.Ic(1):HZ.A11.Ic(2)));  % M/2 cols
act22 = act(ismember(act, HZ.A22.Ic(1):HZ.A22.Ic(2)));  % M/2 cols
M2    = numel(act11);         % = M/2 = 15

perm2 = [act11, act22, nul];
Zv2   = Zw(:, perm2);   % snapshot for plotting

% The cross-level off-diagonal blocks (rank k=5) are now visible here:
% Zv2(1:M2, M2+1:M)  — HZ.A12 rows A11, cols A22  (small, rank ≤ k)
% Zv2(M2+1:M, 1:M2)  — HZ.A21 rows A22, cols A11  (small, rank ≤ k)
fprintf('  Cross-level off-diag [rows 1:%d vs cols %d:%d]: %.2e\n', ...
        M2, M2+1, M, norm(Zv2(1:M2, M2+1:M),'fro'))

%% =========================================================================
%%  PHASE 3 – ROOT / FINAL LQ
%%  Single dense QR on the full working matrix.
%%  Absorbs the residual rank-k cross-level blocks and makes Z·V = [L|0]
%%  exact.  In a full O(N·k²) HSS-ULV this is replaced by a structured
%%  merge step that exploits the rank-k cross blocks.
%% =========================================================================
fprintf('\n─── Phase 3: Root LQ ───────────────────────────────────────────\n')

% Zw' is N×M.  Full QR: Zw' = Qf·Rf  →  Zw·Qf = Rf'
% Rf is N×M;  Rf(1:M,1:M) is upper-tri  →  L = Rf(1:M,1:M)' is lower-tri
% Rf(M+1:N,:) = 0  →  last N-M cols of Rf' are zero  →  Zw·Qf = [L|0]
[Qf, Rf] = qr(Zw');
L        = Rf(1:M, 1:M)';      % M×M lower triangular
Zw       = Zw * Qf;            % = [L | 0] exactly
V        = V  * Qf;

fprintf('  ||triu(L,1)||_F / ||L||_F = %.2e  (lower-tri check)\n', ...
        norm(triu(L,1),'fro') / (norm(L,'fro') + eps))

%% ── Verification ─────────────────────────────────────────────────────────
err = norm(Z*V - [L, zeros(M,N-M)], 'fro') / norm(Z,'fro');
fprintf('\nReconstruction  ||Z·V – [L|0]|| / ||Z||  =  %.3e\n', err)

%% ── Min-norm solution demo ───────────────────────────────────────────────
b     = randn(M,1);
x_ulv = V(:,1:M) * (L \ b);      % U = I  →  U'b = b
x_ref = pinv(Z) * b;

fprintf('\nMin-norm solve:\n')
fprintf('  ||Z·x – b||           =  %.3e\n', norm(Z*x_ulv - b))
fprintf('  ||x_ulv||  =  %.6f   ||x_ref||  =  %.6f\n', norm(x_ulv), norm(x_ref))
fprintf('  ||x_ulv – x_ref||     =  %.3e\n', norm(x_ulv - x_ref))

%% =========================================================================
%%  VISUALISATION
%%  Panel layout:
%%    [Step 0 | Step 1 | Step 2]
%%    [Step 3 | L      | V     ]
%% =========================================================================
data  = {Z,   Zv1,  Zv2,  Zw,  L,  V};
titl  = { ...
    'Step 0 – Z  (original)', ...
    {'Step 1 – after Leaf LQ', 'cols: [active (1:M) | null (M+1:N)]'}, ...
    {'Step 2 – after Level-1 LQ', 'cols: [A11-act | A22-act | null]'}, ...
    {'Step 3 – after Root LQ', 'Z·V = [L | 0]  (exact)'}, ...
    'L  (M×M lower triangular)', ...
    'V  (N×N accumulated right unitary)'};

figure('Name','HSS-ULV Progress','Units','normalized','Position',[.01 .05 .96 .87])
for s = 1:6
    ax = subplot(2,3,s);
    imagesc(abs(data{s})); colorbar; colormap(ax, 'parula'); hold on
    title(titl{s}, 'FontWeight','bold','FontSize',9)
    xlabel('col'); ylabel('row')

    % White dashed lines mark structural boundaries (visualization only)
    switch s
        case 2   % after Leaf LQ
            plot([M+.5 M+.5], [.5 M+.5], 'w--', 'LineWidth', 1.8)   % act|null

        case 3   % after Level-1 LQ
            plot([M2+.5 M2+.5], [.5 M+.5],  'w--', 'LineWidth', 1.8)  % A11|A22 (cols)
            plot([M+.5  M+.5],  [.5 M+.5],  'w--', 'LineWidth', 1.8)  % act|null
            plot([.5 M+.5],     [M2+.5 M2+.5],'w--','LineWidth', 1.8)  % A11|A22 (rows)

        case 4   % after Root LQ
            plot([M+.5 M+.5], [.5 M+.5], 'w--', 'LineWidth', 1.8)   % L|zeros
    end
end
sgtitle('HSS-ULV Factorisation  —  Z · V = [L | 0]', 'FontSize', 12, 'FontWeight','bold')

%% =========================================================================
%%  LOCAL FUNCTIONS
%% =========================================================================

function lv = gather_leaves(node)
%GATHER_LEAVES  Collect HSS leaf nodes in left-to-right order.
    lv = struct('rows',{},'cols',{});
    lv = rec(node, lv);
end

function lv = rec(nd, lv)
%REC  Recursive helper for gather_leaves.
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