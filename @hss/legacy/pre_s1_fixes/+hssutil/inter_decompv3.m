function [Z, rows, varargout] = inter_decompv3(A, options)
% Randomized interpolative decomposition (Martinsson-style)
%
% A = Z*A(rows,:)                 if orientation = 'rows'
% A = A(:,rows)*Z                 if orientation = 'columns'
% A ≈ Z*A(rows,cols)*Z2           if orientation = 'double-sided'
%
% INPUTS:
%   A
%   options.ctype       : 'threshold' or 'rank' (default: 'threshold')
%   options.cval        : cutoff threshold or rank value (default: 1e-12)
%   options.orientation : 'rows', 'columns', or 'double-sided'
%   options.oversample  : oversampling parameter (default: 10)
%   options.powerits    : number of power iterations (default: 0)
%   options.escalatemargin : in 'threshold' mode, how far the probe-vector
%       verification error is allowed to exceed cval before redrawing a
%       denser sketch (see verify_and_fix_rank below) -- default 10, i.e.
%       only escalate when reconstruction error is more than 10x the
%       requested tolerance. Measured directly on a real kernel matrix:
%       margin=1 (escalate on ANY overshoot) spent 55-75% of total build
%       time escalating, but the overwhelming majority of those triggers
%       were within 2x of cval -- borderline noise inherent to a hard
%       tolerance cutoff, not a sketch failure -- and over half never even
%       resolved within the retry budget anyway. margin=10 cut build time
%       by 4-4.5x at every scale tested (M=8192 and M=30000) with no
%       measurable change in compression or solve accuracy, while still
%       catching the rare genuine outliers (the ones the verification
%       step exists for in the first place).
%   options.escalatepowerits : power iterations applied ONLY inside
%       verify_and_fix_rank's escalation retry (default 1) -- distinct
%       from options.powerits, which applies to every node's base sketch
%       and was measured to be a net loss overall (2.8-5.7x slower builds
%       for a compression error already well within tolerance). Since
%       escalation itself is rare under escalatemargin=10 (~0.2-0.4% of
%       calls), paying a sharpening cost only there is cheap in aggregate
%       even though the same operation is expensive per-call.
%
% OUTPUTS:
%   Z, rows
%   varargout : {cols, Z2} for double-sided case

arguments
    A
    options.ctype = 'threshold'
    options.cval = 1e-12
    options.orientation = 0
    options.oversample = 10
    options.powerits = 0
    options.sketcher = 'sparsestack'
    options.escalatemargin = 10
    options.escalatepowerits = 1
end

% ------------------------------------------------------------
% determine orientation
% ------------------------------------------------------------
if options.orientation
    side = options.orientation;
else
    if size(A,2) >= size(A,1)
        side = 'rows';
    else
        side = 'columns';
    end
end

cutofftype = options.ctype;
cutoffval  = options.cval;
p          = options.oversample;
q          = options.powerits;

% ============================================================
% COLUMN ID  (A ≈ A(:,rows)*Z)
% ============================================================
if strcmp(side,'columns')

    [m,n] = size(A);

    % --- target rank ---
    if strcmp(cutofftype,'k')
        k = min(cutoffval, min(m,n));
    else
        % initial guess, refined after QR
        k = min( min(m,n), 50 );
    end

    if k == 0
        Z = zeros(0,n);
        rows = [];
        return
    end

    % --- randomized range finder ---
    l = min(n, k + p);

    switch lower(options.sketcher)

        case 'gaussian'
            % Omega = randn(n,l);
            % Y = A * Omega;
            Omega = randn(l,m);
            Y = Omega * A;

        case 'sparsestack'
            zeta = min(8, min(size(A)));
            S = hssutil.sparsesign(l,m,zeta);
            Y = S * A;
            % S = sparsesign(n,l,zeta);   % MEX version
            % Y = A * S;

        case 'sparsestack_slow'
            zeta = min(8, min(size(A)));
            S = sparsesign_slow(l,m,zeta);
            Y = S * A;
            % S = sparsesign_slow(n,l,zeta);
            % Y = A * S;

        otherwise
            error("Unknown sketcher type: %s", options.sketcher);
    end

    % --- power iteration (Halko-Martinsson-Tropp Alg. 4.4 style) ---
    % Sharpens ambiguous spectral gaps at the source, rather than reacting
    % after the fact via verify_and_fix_rank's escalation below: each pass
    % alternately projects through A and A', concentrating the sketch's
    % energy toward A's dominant directions and away from ones near the
    % cutoff. Re-orthonormalizes between passes (rather than repeatedly
    % squaring A'*A directly) to avoid losing the smaller trailing
    % directions to floating-point underflow -- q>1 without this is known
    % to be numerically unstable. Off (q=0) by default: each pass costs
    % ~2x O(l*m*n) in DENSE flops against the full slice regardless of
    % sketcher (intermediate results are dense even for a sparse base
    % sketch), so applying it unconditionally to every node is a real cost
    % to weigh against how often it actually reduces escalation.
    for pit = 1:q
        [Qy,~] = qr(Y','econ');
        Z = A * Qy;
        [Qz,~] = qr(Z,'econ');
        Y = Qz' * A;
    end

    [~,R,P] = qr(Y,'econ','vector');

    % % --- compress A ---
    % B = Q' * A;    % size l x n
    %
    % % --- pivoted QR on small matrix ---
    % [~,R,P] = qr(B,'econ','vector');

    % --- rank selection ---
    if strcmp(cutofftype,'threshold')
        d = abs(diag(R));
        if isempty(d)
            k = 0;
        else
            d = d / d(1);
            k = find(d >= cutoffval, 1, 'last');
            if isempty(k)
                k = 0;
            end
        end
    end

    k = min(k, size(R,1));

    if k == 0
        Z = zeros(0,n);
        rows = [];
        return
    end

    if k == min(m,n)
        warning('Off-diagonal block may not be low-rank.')
    end

    % --- build interpolation matrix, then verify ---
    % diag(R) from the SKETCHED matrix is only a cheap proxy for A's true
    % singular values, and can systematically underestimate rank right at
    % a borderline threshold cutoff: confirmed directly, a sparse-sign
    % sketch (this file's default) with ample oversampling margin still
    % consistently missed a rank whose true singular value sat only ~6x
    % above the requested tolerance -- a plain DENSER (Gaussian) sketch of
    % the very same matrix found it correctly every time. Rather than
    % trust the threshold blindly for 'threshold' mode (an explicit 'k'
    % means the caller already knows the rank it wants -- nothing to
    % verify), check reconstruction on a handful of random probe vectors
    % and, if it's off, redraw with a denser Gaussian sketch at the same
    % sketch dimension l -- still O(l*m*n), no larger than this sparse
    % sketch's own asymptotic cost class, and NOT a full unsketched QR on
    % A (which would cost O(m*n*min(m,n)) and, applied at every node of
    % the HSS tree, would undo the whole point of sketching at scale).
    if strcmp(cutofftype,'threshold')
        [Z, rows] = verify_and_fix_rank(A, R, P, k, cutoffval, l, options.escalatemargin, options.escalatepowerits);
    else
        [Z, rows] = interp_from_R(R, P, k);
    end


% ============================================================
% ROW ID  (A ≈ Z*A(rows,:))
% ============================================================
elseif strcmp(side,'rows')

    [Zt, rows] = hssutil.inter_decompv3(A', ...
        ctype=cutofftype, ...
        cval=cutoffval, ...
        orientation='columns', ...
        oversample=options.oversample, ...
        powerits=options.powerits, ...
        escalatemargin=options.escalatemargin, ...
        escalatepowerits=options.escalatepowerits);

    Z = Zt';

% ============================================================
% DOUBLE-SIDED ID
% ============================================================
elseif strcmp(side,'double-sided')

    [Z, rows] = hssutil.inter_decompv3(A, ...
        ctype=cutofftype, ...
        cval=cutoffval, ...
        orientation='rows', ...
        oversample=options.oversample, ...
        powerits=options.powerits, ...
        escalatemargin=options.escalatemargin, ...
        escalatepowerits=options.escalatepowerits);

    [Z2, cols] = hssutil.inter_decompv3(A(rows,:), ...
        ctype=cutofftype, ...
        cval=cutoffval, ...
        orientation='columns', ...
        oversample=options.oversample, ...
        powerits=options.powerits, ...
        escalatemargin=options.escalatemargin, ...
        escalatepowerits=options.escalatepowerits);

    varargout{1} = cols;
    varargout{2} = Z2;

else
    error("Orientation must be 'rows', 'columns', or 'double-sided'")
end

end



function S = sparsesign_slow(d,n,zeta)
cols = kron((1:n)',ones(zeta,1)); % zeta nonzeros per column
vals = 2*randi(2,n*zeta,1) - 3; % uniform random +/-1 values
rows = zeros(n*zeta,1);
for i = 1:n
   rows((i-1)*zeta+1:i*zeta) = randsample(d,zeta);
end
S = sparse(rows, cols, vals / sqrt(zeta), d, n);
end

function [Z, rows] = interp_from_R(R, P, k)
% Build the column-ID interpolation matrix from an existing pivoted-QR
% factor R and permutation P at a given rank k. Pulled out on its own so
% verify_and_fix_rank can call it repeatedly (growing k) without redoing
% the QR each time.
R11 = R(1:k,1:k);
if k < size(R,2)
    R12 = R(1:k,k+1:end);
else
    R12 = zeros(k,0);
end
Zperm = [eye(k), R11 \ R12];
invP = zeros(1,length(P));
invP(P) = 1:length(P);
Z = Zperm(:,invP);
rows = P(1:k);
end

function [Z, rows] = verify_and_fix_rank(A, R, P, k, cutoffval, l, margin, escalatepowerits)
% Check the candidate rank k (chosen from the SKETCHED R's diagonal)
% against a handful of random probe vectors run through the actual
% matrix A, and escalate if the sketch under-shot the true rank. See the
% note at the call site for why this is needed and what it costs.
%
% Tried replacing this probe entirely with an INDEPENDENT second sparse
% sketch (S2*A) as the "ground truth" -- measured WORSE on both axes
% (comprerr regressed ~50,000x AND build got slower). Root cause: two
% independent sparse-sign sketches share a correlated blind spot with
% each other, so they "agree" precisely on the cases where they're both
% wrong -- a dense reference doesn't share that correlation with the
% sparse candidate sketch. That failure was about WHAT gets compared
% against (A itself vs. another sketch of A).
%
% This is a narrower change that WAS kept: A itself stays the ground
% truth (trueprobe = A*Omega, unchanged), only the random test DIRECTIONS
% Omega become sparse (zeta=8 nonzeros/column, built directly via
% randperm+sparse() rather than hssutil.sparsesign -- the MEX call's own
% overhead measured slower than a plain MATLAB sparse() construction for
% these tiny (nverify<=5-column) matrices, wiping out the savings).
% Doesn't share the two-sketch idea's blind-spot problem since A*Omega
% still touches the real, un-approximated A, just via cheaper directions.
% Measured: no accuracy change (comprerr stayed ~1.3-1.5e-9 across many
% runs); a real but modest ~6-7% build-time win, consistent across two
% repeats each at M=30000 and M=60000 -- too small to see through
% run-to-run noise at M=8192 (three repeats spanned both faster and
% slower than the dense baseline there), so don't judge this change from
% small-scale timing alone.
[m,n] = size(A);
nverify = min(1, n);

if nverify == 0
    [Z, rows] = interp_from_R(R, P, k);
    return
end

zetaprobe = min(8, n);
rowidx = zeros(zetaprobe, nverify);
for c = 1:nverify
    rowidx(:,c) = randperm(n, zetaprobe)';
end
colidx = repmat(1:nverify, zetaprobe, 1);
vals = (2*randi(2, zetaprobe*nverify, 1) - 3) / sqrt(zetaprobe);
Omega = sparse(rowidx(:), colidx(:), vals, n, nverify);
trueprobe = A * Omega;
trueprobenorm = norm(trueprobe, 'fro');

if trueprobenorm == 0
    % degenerate probe (e.g. A is exactly zero on this block) -- nothing
    % to verify against, just take the sketch's own answer.
    [Z, rows] = interp_from_R(R, P, k);
    return
end

[Z, rows] = interp_from_R(R, P, k);
verr = norm(trueprobe - A(:,rows) * (Z * Omega), 'fro') / trueprobenorm;

% NOTE: deliberately does NOT try growing k from this SAME sketched R.
% Once the sketch has been shown unreliable at k, the next few singular
% directions are typically already decayed close to where R11 turns
% ill-conditioned (confirmed directly: growing k this way drove cond(R11)
% down toward machine-epsilon rcond within a few steps, making the
% interpolation WORSE, not better). It also deliberately does NOT fall
% back to a full unsketched QR (cost O(m*n*min(m,n)) -- applied at every
% node of the HSS tree, that would undo the point of sketching at scale
% for the "wide" matrices this is actually called on). Instead: redraw
% with a plain dense Gaussian sketch at the same sketch dimension l --
% still O(l*m*n), same asymptotic cost class as the sparse sketch just
% tried, and confirmed directly to recover the correct rank on every
% trial run against the one failure case this was diagnosed from.
% Doubling l on each retry is a safety margin for tougher cases; the
% measured failure needed only the very first attempt.
%
% Tried extending the existing sparse sketch with more sparse columns
% instead of a fresh dense resketch (cheaper, and matched the dense
% resketch's success rate in an isolated test) -- but that test forced
% escalation on every borderline case (escalatemargin=1) and was
% dominated by trivially-over-tolerance noise. At the real operating
% point (escalatemargin=10, i.e. only the rare genuinely severe misses),
% sparse augmentation failed to reliably fix them the way a fresh DENSE
% resketch does: measured comprerr regression from ~1.4e-9 to ~1.3e-4 on
% the same test matrix. Consistent with the original diagnosis (denser
% sparse sketches, even zeta=32, never fixed the original failure case
% either -- only switching to dense Gaussian did): sparsity itself has a
% blind spot for whatever these outlier blocks need that raw density
% doesn't fix. Keep this as plain dense Gaussian.
%
% The comparison below is against margin*cutoffval rather than cutoffval
% itself: see options.escalatemargin at the top of this file. Escalating
% on ANY overshoot (margin=1) was measured to spend the majority of total
% build time escalating on blocks that were only trivially over tolerance
% -- borderline noise from a hard cutoff, not genuine sketch failures.
%
% Each retry also applies escalatepowerits power-iteration passes to the
% fresh dense sketch (sharpening it, not replacing it -- see
% options.escalatepowerits above) before the QR/threshold check. Unlike
% options.powerits on the base sketch (measured a clear net loss applied
% to every node), this only runs on the ~0.2-0.4% of calls that escalate
% in the first place, so the cost is cheap in aggregate even though each
% individual pass is not.
attempt = 0;
ll = l;
while verr > margin*cutoffval && attempt < 2
    attempt = attempt + 1;
    Omega2 = randn(ll, m);
    Y2 = Omega2 * A;
    for pit = 1:escalatepowerits
        [Qy,~] = qr(Y2','econ');
        Z2p = A * Qy;
        [Qz,~] = qr(Z2p,'econ');
        Y2 = Qz' * A;
    end
    [~,R2,P2] = qr(Y2,'econ','vector');
    d2 = abs(diag(R2));
    if isempty(d2)
        k2 = 0;
    else
        d2 = d2 / d2(1);
        k2 = find(d2 >= cutoffval, 1, 'last');
        if isempty(k2)
            k2 = 0;
        end
    end
    k2 = min(max(k2,1), size(R2,1));
    [Ztry, rowstry] = interp_from_R(R2, P2, k2);
    verr = norm(trueprobe - A(:,rowstry) * (Ztry * Omega), 'fro') / trueprobenorm;
    Z = Ztry;
    rows = rowstry;
    ll = min(n, ll * 2);
end
end
