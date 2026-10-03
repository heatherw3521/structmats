function [x,H] = hss_ulvminnormsolve(H,b)
% ULV solve with a wide or square multi-level HSS matrix:
%     x = argmin ||x||_2  subject to  H*x = b
% (for square nonsingular H, the unique solution). H\b calls this for both
% shapes; the factors are kept in H's factor cache (see @hss/hss.m), so
% users never call this directly.
%
%   x = hss_ulvminnormsolve(H, b)    factor, then solve (b may have several columns)
%   F = hss_ulvminnormsolve(H)       factor only; F is a struct of stored factors
%   x = hss_ulvminnormsolve(F, b)    solve with stored factors, O(n r) per solve
%   [x, H] = hss_ulvminnormsolve(H, b) returns H unchanged (signature used by mldivide)
%
% Public entry points: H\b (mldivide.m), minnorm.m, tikhonov.m, which cache
% the struct F returned by the factor-only call.
%
% Requirements:
%   * H has full row rank. 'hss_ulvminnormsolve:rankDeficient' is raised when
%     the sizes show that it fails (a leaf with more decoupled rows than
%     columns, or a tall root block). A deficiency the sizes do not show (a
%     singular L_tau, or a rank-deficient wide root block) is not detected.
%   * no row basis is wider than the block it spans: k_tau <= m_tau at every
%     (current) leaf ('hss_ulvminnormsolve:improperRanks', always checked).
%     Nested-ID constructions (hss_constructor) always satisfy this.
%
% Changes from @hss/legacy/hss_ulvminnormsolve.m (memo 1, HSS_Min_Norm_Solve_Memo_v2.pdf):
%   1. Root solve applies the SVD factors of the root block to b. The legacy
%      pinv(H.D)*b forms the pseudo-inverse explicitly, which is not backward
%      stable: its residual grows like cond(root)*eps*norm(b) (memo 1, Sec. 6.4).
%   2. The checks above replace failures inside an index.
%   3. Size-reduced width p' = min(n_tau, l_tau + m_tau) instead of l_tau + m_tau.
%      The legacy fixed width needs l_tau + m_tau <= n_tau ("slack"); without it,
%      YD(:,1:(l+m)) indexes past the n_tau columns the LQ factor has, which is
%      why the legacy code checks the slack and merges tree levels (levelup)
%      until it holds. With the relaxed width nothing is discarded when the
%      slack is short, and the only size condition left, t_tau = m_tau - k_tau
%      <= p'_tau (L_tau square), holds whenever H has full row rank (memo 1,
%      Lemma 4.4 and Prop. 4.8). No level is ever merged; square H works too.
%   4. Economic QR in the size reduction: Omega_tau is p' x n_tau instead of
%      n_tau x n_tau, so time and memory drop by about n_tau/(l_tau + m_tau)
%      for very wide matrices. (The QL of U and the decoupling LQ stay full:
%      P_tau and Q_tau must be square.)
%   5. Factor/solve split: every level's factors (Omega, P, Q, L, the
%      transformed matrix used for the RHS update) are stored, so a new
%      right-hand side costs triangular solves, one HSS product per level and
%      small unitary products, without refactoring. Several right-hand sides
%      (columns of b) are solved together.
%   6. A leaf with l_tau + m_tau >= n_tau (always the case for square H) skips
%      the size reduction: it would discard nothing, so Omega_tau = I. This is
%      what makes the same code the square solver too (it replaced
%      hss_ulvvecsolve.m, now in @hss/legacy/, as the square path of H\b).
%
% Notation in this file follows the code: a leaf has m rows and n columns,
% k = rank of its row basis (Z), l = rank of its column basis (Y); in the
% memo these are n_tau, p_tau, k_tau, l_tau.

if nargin == 1
    x = factor_level(H);
    return
end
if isstruct(H)
    x = solve_level(H, b);
    return
end
x = solve_level(factor_level(H), b);
end


% ======================================================================
% factorization
% ======================================================================
function F = factor_level(H)
% Factor one recursion level and, recursively, the reduced matrix.
H.factorcache = [];                % work on a cache-free copy (no set-method churn)
if H.isleaf
    F = root_factor(H.D);
    return
end
acc = struct('m', [], 'n', [], 'r', [], 't', [], 'Om', {{}}, 'P', {{}}, ...
             'Q', {{}}, 'Lt', {{}});
[HULV, acc] = down_recurse(H, acc);

nl = numel(acc.m);
ro = [0, cumsum(acc.m)];          % row offsets of the leaves (this level)
co = [0, cumsum(acc.r)];          % column offsets in the size-reduced coordinates
xo = [0, cumsum(acc.n)];          % column offsets in the original coordinates
forced_rows = false(ro(end), 1);
forced_cols = false(co(end), 1);
for i = 1:nl
    forced_rows(ro(i) + (1:acc.t(i))) = true;
    forced_cols(co(i) + (1:acc.t(i))) = true;
end

F.isroot = false;
F.nleaves = nl;
F.m = acc.m; F.n = acc.n; F.r = acc.r; F.t = acc.t;
F.ro = ro; F.co = co; F.xo = xo;
F.Om = acc.Om; F.P = acc.P; F.Q = acc.Q; F.Lt = acc.Lt;
F.keep_rows = find(~forced_rows);
F.free_cols = find(~forced_cols);
F.HULV = HULV;                    % transformed matrix of this level (RHS update)

% reduced system: drop decoupled rows / forced columns, merge sibling leaves
F.next = factor_level(levelup(discard(HULV), 1, 1));
end


function F = root_factor(D)
% Root block: a single dense block, square or wide, of full row rank.
[m, n] = size(D);
F.isroot = true;
F.size = [m, n];
if m > n
    error('hss_ulvminnormsolve:rankDeficient', ...
        ['the reduced root block is %d x %d (more rows than columns), so H ' ...
         'does not have full row rank.'], m, n);
end
if m == n
    [F.L, F.U, F.p] = lu(D, 'vector');
    F.kind = 'lu';
else
    % minimum-norm solve from the SVD factors, with pinv's default rank
    % tolerance; never form pinv(D) itself (not backward stable)
    [U, S, V] = svd(D, 'econ');
    s = diag(S);
    if isempty(s) || s(1) == 0
        r = 0;
    else
        r = sum(s > max(m, n)*eps(s(1)));
    end
    F.U = U(:, 1:r); F.s = s(1:r); F.V = V(:, 1:r);
    % left singular vectors of the dropped directions: empty when the root
    % has full row rank; otherwise root_solve uses them to tell a
    % consistent b from one with no solution (audit B11)
    F.Uperp = U(:, r+1:m);
    F.kind = 'svd';
end
end


function [H, acc] = down_recurse(H, acc)
% Size reduction and decoupling of every leaf pair below H. Collects the
% per-leaf factors in left-to-right leaf order and returns the transformed
% matrix (sizes and column index ranges refreshed to the reduced widths).
if H.A11.isleaf
    % leaf 1 = A11: rows m1, cols n1, row basis H.A12.Z (k1), col basis H.A21.Y (l1)
    % leaf 2 = A22: rows m2, cols n2, row basis H.A21.Z (k2), col basis H.A12.Y (l2)
    m1 = H.A11.sz(1);  n1 = H.A11.sz(2);
    m2 = H.A22.sz(1);  n2 = H.A22.sz(2);
    k1 = size(H.A12.Z, 2);  l1 = size(H.A21.Y, 1);
    k2 = size(H.A21.Z, 2);  l2 = size(H.A12.Y, 1);
    check_ranks(H.A11.Ir, k1, m1);
    check_ranks(H.A22.Ir, k2, m2);

    % ---- size reduction (change 3: width min(n, l+m); change 4: economic QR)
    % [V^*; D] = L*Om, Om with orthonormal rows; keep all r = min(n, l+m)
    % columns of L.  Columns outside Om's row space are zero in the whole
    % block column (memo 1, Lemma 2.1) and are dropped.  Change 6: if
    % l + m >= n nothing would be dropped, so skip it (Om = [] means I).
    if l1 + m1 < n1
        [YD1, Om1] = lq_econ([H.A21.Y; H.A11.D]);
        H.A21.Y = YD1(1:l1, :);   H.A11.D = YD1(l1+1:end, :);
    else
        Om1 = [];
    end
    if l2 + m2 < n2
        [YD2, Om2] = lq_econ([H.A12.Y; H.A22.D]);
        H.A12.Y = YD2(1:l2, :);   H.A22.D = YD2(l2+1:end, :);
    else
        Om2 = [];
    end
    r1 = size(H.A11.D, 2);    r2 = size(H.A22.D, 2);

    % ---- decoupling, step 1: QL of the row bases, Z = P*[0; Zhat]
    [P1, Zhat1] = ql(H.A12.Z);
    [P2, Zhat2] = ql(H.A21.Z);
    H.A12.Z = Zhat1;
    H.A21.Z = Zhat2;
    H.A11.D = P1'*H.A11.D;
    H.A22.D = P2'*H.A22.D;

    % ---- decoupling, step 2: LQ of the top t = m - k rows (they now couple to
    % nothing outside the leaf); apply Q to the leaf's whole block column
    t1 = m1 - k1;  t2 = m2 - k2;
    check_width(H.A11.Ir, t1, r1);
    check_width(H.A22.Ir, t2, r2);
    [Lt1, Q1] = lq(H.A11.D(1:t1, :));
    [Lt2, Q2] = lq(H.A22.D(1:t2, :));
    H.A11.D(1:t1, :) = Lt1;
    H.A22.D(1:t2, :) = Lt2;
    H.A11.D(t1+1:end, :) = H.A11.D(t1+1:end, :)*Q1';
    H.A22.D(t2+1:end, :) = H.A22.D(t2+1:end, :)*Q2';
    H.A21.Y = H.A21.Y*Q1';
    H.A12.Y = H.A12.Y*Q2';

    % ---- store the factors of both leaves (leaf order: 1, 2)
    acc.m(end+1:end+2) = [m1, m2];
    acc.n(end+1:end+2) = [n1, n2];
    acc.r(end+1:end+2) = [r1, r2];
    acc.t(end+1:end+2) = [t1, t2];
    acc.Om(end+1:end+2) = {Om1, Om2};
    acc.P(end+1:end+2)  = {P1, P2};
    acc.Q(end+1:end+2)  = {Q1, Q2};
    acc.Lt(end+1:end+2) = {Lt1(:, 1:t1), Lt2(:, 1:t2)};

    % ---- refresh size/Ic bookkeeping to the reduced widths (needed by
    % hss_matvec, and by the ancestors' index shifting)
    H.A11.sz = size(H.A11.D);
    H.A22.sz = size(H.A22.D);
    H.A12.sz = [size(H.A12.Z,1), size(H.A12.Y,2)];
    H.A21.sz = [size(H.A21.Z,1), size(H.A21.Y,2)];
    H.A11.Ic = [H.Ic(1), H.Ic(1)+H.A11.sz(2)-1];
    H.A22.Ic = [H.A11.Ic(2)+1, H.A11.Ic(2)+H.A22.sz(2)];
    H.A12.Ic = H.A22.Ic;  H.A12.Ir = H.A11.Ir;
    H.A21.Ic = H.A11.Ic;  H.A21.Ir = H.A22.Ir;
    H.Ic = [H.A11.Ic(1), H.A22.Ic(2)];
    H.sz = [H.A11.sz(1)+H.A22.sz(1), H.A11.sz(2)+H.A22.sz(2)];
else
    [H.A11, acc] = down_recurse(H.A11, acc);
    [H.A22, acc] = down_recurse(H.A22, acc);
    H.sz = [H.A11.sz(1)+H.A22.sz(1), H.A11.sz(2)+H.A22.sz(2)];
    H.Ic = [H.A11.Ic(1), H.A22.Ic(2)];
    % off-diagonal Ic/size at this level: A12 spans rows of A11, cols of A22
    H.A12.Ic = H.A22.Ic;  H.A12.Ir = H.A11.Ir;
    H.A21.Ic = H.A11.Ic;  H.A21.Ir = H.A22.Ir;
    H.A12.sz = [H.A11.sz(1), H.A22.sz(2)];
    H.A21.sz = [H.A22.sz(1), H.A11.sz(2)];
    % H.Ir is unchanged
end
end


function check_ranks(Ir, k, m)
% change 2: a row basis wider than its leaf would make t = m - k negative
if k > m
    error('hss_ulvminnormsolve:improperRanks', ...
        ['the leaf with rows %s has a row basis of rank %d but only %d rows. ' ...
         'Make the representation proper (replace the basis by an orthonormal ' ...
         'basis of its range and absorb the rest into the parent translation ' ...
         'and coupling matrices); hss_constructor never produces this.'], ...
        mat2str(Ir), k, m);
end
end


function check_width(Ir, t, r)
% change 2: the t decoupled rows of a leaf live on its r size-reduced columns
% only; if t > r they are linearly dependent, so H lacks full row rank
if t > r
    error('hss_ulvminnormsolve:rankDeficient', ...
        ['the leaf with rows %s has %d rows that couple to no other leaf but ' ...
         'only %d columns, so H does not have full row rank.'], mat2str(Ir), t, r);
end
end


% ======================================================================
% solve with stored factors
% ======================================================================
function X = solve_level(F, B)
if F.isroot
    X = root_solve(F, B);
    return
end
if size(B, 1) ~= F.ro(end)
    error('hss_ulvminnormsolve:dimension', ...
        'b has %d rows; the matrix has %d.', size(B, 1), F.ro(end));
end
ns = size(B, 2);
W = zeros(F.co(end), ns);          % solution in the size-reduced, decoupled coordinates
for i = 1:F.nleaves
    rows = F.ro(i) + (1:F.m(i));
    B(rows, :) = F.P{i}'*B(rows, :);                 % row compression of b
    t = F.t(i);
    if t > 0                                         % forced solve, L_tau lower triangular
        W(F.co(i) + (1:t), :) = F.Lt{i} \ B(F.ro(i) + (1:t), :);
    end
end
% reduced right-hand side: b_hat - H_hat*[X1; 0] on the kept rows (all leaves'
% forced unknowns contribute, through the nested bases -- one HSS product)
Bbar = B - F.HULV*W;
W(F.free_cols, :) = solve_level(F.next, Bbar(F.keep_rows, :));
% reconstruction: x_tau = Om_tau' * Q_tau' * [X1; X2]
X = zeros(F.xo(end), ns);
for i = 1:F.nleaves
    Xi = F.Q{i}'*W(F.co(i) + (1:F.r(i)), :);
    if ~isempty(F.Om{i}), Xi = F.Om{i}'*Xi; end      % Om = [] : no size reduction
    X(F.xo(i) + (1:F.n(i)), :) = Xi;
end
end


function X = root_solve(F, B)
if size(B, 1) ~= F.size(1)
    error('hss_ulvminnormsolve:dimension', ...
        'b has %d rows; the matrix has %d.', size(B, 1), F.size(1));
end
if strcmp(F.kind, 'lu')
    X = F.U \ (F.L \ B(F.p, :));
else
    if ~isempty(F.Uperp)
        check_consistent(F, B);
    end
    X = F.V*((F.U'*B)./F.s);
end
end


function check_consistent(F, B)
% The root block has numerical rank r < m, so H does not have full row rank
% (a dependence between rows of different leaves reaches the root). The
% part of the reduced right-hand side in the dropped directions decides:
%  * about zero: H*x = b has solutions, and the SVD solve returns the
%    minimum-norm one exactly (all earlier steps are exact reformulations
%    of a solvable problem) -> warning, return x;
%  * not zero: H*x = b has no solution. The forced solves of the earlier
%    levels have already imposed their rows exactly, so the result would be
%    neither a solution nor the least-squares minimum-norm solution
%    -> error.
tol = sqrt(eps);
nb = sqrt(sum(abs(B).^2, 1));
nd = sqrt(sum(abs(F.Uperp'*B).^2, 1));
ratio = max(nd ./ max(nb, realmin));
m = F.size(1); r = m - size(F.Uperp, 2);
if ratio > tol
    error('hss_ulvminnormsolve:inconsistent', ...
        ['H does not have full row rank (the reduced root block has rank %d of %d) and ' ...
         'H*x = b has no solution (relative inconsistency %.1e). The minimum-norm solver ' ...
         'needs a consistent system; for a regularized least-squares solution use ' ...
         'tikhonov(H, b, lambda).'], r, m, ratio);
end
warning('hss_ulvminnormsolve:rankDeficientRoot', ...
    ['H does not have full row rank (the reduced root block has rank %d of %d); b is ' ...
     'consistent (relative inconsistency %.1e), so the result is the minimum-norm ' ...
     'solution, but it is sensitive to perturbations of b.'], r, m, ratio);
end


% ======================================================================
% reduced matrix: discard decoupled rows / forced columns, merge leaves
% (unchanged from the legacy file)
% ======================================================================
function H = discard(H)
if H.A11.isleaf
    [m1,n1] = size(H.A11.D);
    [m2,n2] = size(H.A22.D);
    k1 = size(H.A12.Z,2);
    k2 = size(H.A21.Z,2);

    H.A12.Z = H.A12.Z(end-k1+1:end,:);
    H.A12.Y = H.A12.Y(:,m2-k2+1:end);
    H.A12.sz = [size(H.A12.Z,1), size(H.A12.Y,2)];
    H.A12.isleaf = 1;

    H.A21.Z = H.A21.Z(end-k2+1:end,:);
    H.A21.Y = H.A21.Y(:,m1-k1+1:end);
    H.A21.sz = [size(H.A21.Z,1), size(H.A21.Y,2)];
    H.A21.isleaf = 1;

    H.A11.D = H.A11.D(m1-k1+1:end, m1-k1+1:end);
    H.A11.sz = size(H.A11.D);
    H.A22.D = H.A22.D(m2-k2+1:end, m2-k2+1:end);
    H.A22.sz = size(H.A22.D);
else
    H.A11 = discard(H.A11);
    H.A22 = discard(H.A22);
end
end


function [H,rstart,cstart] = levelup(H,rstart,cstart)
if H.A11.isleaf
    H.D = [H.A11.D, H.A12.Z*H.A12.lrcomponent*H.A12.Y; H.A21.Z*H.A21.lrcomponent*H.A21.Y, H.A22.D];
    H.A11 = [];
    H.A22 = [];
    H.A12 = [];
    H.A21 = [];
    H.isleaf = 1;
    H.sz = size(H.D);
    H.levelcount = H.levelcount-1;
elseif H.A11.A11.isleaf
    [H,rstart,cstart] = onerebuild(H.A12,H,rstart,cstart);
    [H,rstart,cstart] = onerebuild(H.A21,H,rstart,cstart);
    [H,rstart,cstart] = onerebuild(H.A11,H,rstart,cstart);
    [H,rstart,cstart] = onerebuild(H.A22,H,rstart,cstart);
    H.A12.Ir = H.A11.Ir;
    H.A12.Ic = H.A22.Ic;
    H.A21.Ir = H.A22.Ir;
    H.A21.Ic = H.A11.Ic;
    H.sz = H.A11.sz + H.A22.sz;
    H.Ir = [H.A11.Ir(1) H.A22.Ir(2)];
    H.Ic = [H.A11.Ic(1) H.A22.Ic(2)];
    H.levelcount = H.levelcount-1;
else
    [H.A11,rstart,cstart] = levelup(H.A11,rstart,cstart);
    [H.A22,rstart,cstart] = levelup(H.A22,rstart,cstart);
    H.A12.Ir = H.A11.Ir;
    H.A12.Ic = H.A22.Ic;
    H.A21.Ir = H.A22.Ir;
    H.A21.Ic = H.A11.Ic;
    H.sz = H.A11.sz + H.A22.sz;
    H.Ir = [H.A11.Ir(1) H.A22.Ir(2)];
    H.Ic = [H.A11.Ic(1) H.A22.Ic(2)];
    H.levelcount = H.levelcount-1;
end
end


function [Hparent,rstart,cstart] = onerebuild(H,Hparent,rstart,cstart)
% UNCHANGED from hss_ulvvecsolve.m.
if H.isdiag
    H.D = [H.A11.D, H.A12.Z*H.A12.lrcomponent*H.A12.Y; H.A21.Z*H.A21.lrcomponent*H.A21.Y,H.A22.D];
    H.A11 = [];
    H.A12 = [];
    H.A21 = [];
    H.A22 = [];
    H.sz = size(H.D);
    H.levelcount = H.levelcount-1;
    H.isleaf = 1;
    H.Ir = [rstart, rstart + H.sz(1)-1];
    H.Ic = [cstart, cstart + H.sz(2)-1];
    rstart = rstart + H.sz(1);
    cstart = cstart + H.sz(2);
    if mod(H.rowtreeindex,2) ==1
        Hparent.A11 = H;
    else
        Hparent.A22 = H;
    end
elseif mod(H.rowtreeindex,2) ==1
    Z1 = Hparent.A11.A12.Z;
    Z2 = Hparent.A11.A21.Z;
    Hparent.A12.Z = blkdiag(Z1,Z2)*Hparent.A12.Z;

    Y1 = Hparent.A22.A21.Y;
    Y2 = Hparent.A22.A12.Y;
    Hparent.A12.Y = Hparent.A12.Y*blkdiag(Y1,Y2);
    H.levelcount = H.levelcount-1;
    Hparent.A12.sz = [size(Hparent.A12.Z,1),size(Hparent.A12.Y,2)];

elseif mod(H.rowtreeindex,2) == 0
    Z1 = Hparent.A22.A12.Z;
    Z2 = Hparent.A22.A21.Z;
    Hparent.A21.Z = blkdiag(Z1,Z2)*Hparent.A21.Z;

    Y1 = Hparent.A11.A21.Y;
    Y2 = Hparent.A11.A12.Y;
    Hparent.A21.Y = Hparent.A21.Y*blkdiag(Y1,Y2);
    H.levelcount = H.levelcount-1;
    Hparent.A21.sz = [size(Hparent.A21.Z,1),size(Hparent.A21.Y,2)];
end

end


% ======================================================================
% dense kernels.  Convention: lq returns A = L*Q (Q with orthonormal rows),
% so the stored factors are the conjugate transposes of the memo's
% Omega_tau, Q_tau, and the reconstruction applies Om' and Q'.
% ======================================================================
function [Q, L] = ql(A)
% A = Q*L, Q square unitary, L zero except its last size(A,2) rows
[Q_r, R] = qr(flipud(fliplr(A)));
L = flipud(fliplr(R));
Q = flipud(fliplr(Q_r));
end

function [L, Q] = lq(A)
% full: Q square (needed for the decoupling, whose free columns are its complement)
[Q_t, R] = qr(A');
L = R';
Q = Q_t';
end

function [L, Q] = lq_econ(A)
% economic: A (a x n) = L*Q, Q of size min(a,n) x n with orthonormal rows
[Q_t, R] = qr(A', 0);
L = R';
Q = Q_t';
end
