function [Z, rows, varargout] = inter_decompv3(A, options)
%INTER_DECOMPV3  Randomized interpolative decomposition (ID).
%
% A = Z*A(rows,:)                 if orientation = 'rows'
% A = A(:,rows)*Z                 if orientation = 'columns'
% A ≈ Z*A(rows,cols)*Z2           if orientation = 'double-sided'
%
% INPUTS:
%   A
%   options.ctype       : 'threshold' (relative tolerance) or 'k' (fixed rank);
%                         default 'threshold'
%   options.cval        : the tolerance or the rank (default: 1e-12)
%   options.orientation : 'rows', 'columns', or 'double-sided'
%   options.oversample  : oversampling parameter p (default: 10)
%   options.powerits    : power iterations on every sketch (default: 0)
%   options.sketcher    : 'sparsestack' (sparse sign, default), 'gaussian'
%   options.escalatemargin : in 'threshold' mode, the ID is checked on a
%       random probe vector; if its relative error exceeds escalatemargin*cval
%       the block is resketched with a dense Gaussian sketch (default 10, so
%       that errors just above the tolerance do not trigger a redraw)
%   options.escalatepowerits : power iterations applied only to those dense
%       resketches (default 1)
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
if ~(strcmp(cutofftype, 'threshold') || strcmp(cutofftype, 'k'))
    error('hssutil:inter_decompv3:ctype', 'ctype must be ''threshold'' or ''k'' (got ''%s'').', char(cutofftype));
end
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
    [R, P] = sketch_qr(A, l, options.sketcher, q);
    if strcmp(cutofftype,'threshold')
        k = rank_from_R(R, cutoffval);
        % Adaptive sketch size: a sketch with l rows cannot detect a rank
        % above l. While the detected rank leaves fewer than p rows of
        % oversampling, double l and resketch; at l = min(m,n) the sketch
        % carries the whole row space, so the loop always ends.
        while k > l - p && l < min(m,n)
            l = min(min(m,n), 2*l);
            [R, P] = sketch_qr(A, l, options.sketcher, q);
            k = rank_from_R(R, cutoffval);
        end
    end

    k = min(k, size(R,1));

    if k == 0
        Z = zeros(0,n);
        rows = [];
        return
    end

    if k == min(m,n)
        warning('hssutil:notLowRank', 'Off-diagonal block may not be low-rank (numerical rank = min(m,n)).')
    end

    % --- build interpolation matrix, then verify ---
    % diag(R) of the sketched matrix only estimates the singular values of
    % A, and a sparse sketch can underestimate the rank near the threshold.
    % In 'threshold' mode the ID is therefore checked against A on a random
    % probe vector and, if it is off, redrawn with a dense Gaussian sketch of
    % the same size (still O(l m n), unlike an unsketched pivoted QR). A
    % fixed rank ('k') needs no check.
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
% Column-ID interpolation matrix from a pivoted-QR factor R and permutation P
% at rank k.
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
% Check the candidate rank k (chosen from the sketched R's diagonal) against
% A itself on a random probe vector Omega (sparse, zeta = 8 nonzeros per
% column), and redraw with a dense Gaussian sketch if the sketch missed part
% of the range. The reference A*Omega uses A, not a second sketch, so the
% check does not share the sketch's blind spots.
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
    % A*Omega = 0 (e.g. a zero block): nothing to check against
    [Z, rows] = interp_from_R(R, P, k);
    return
end

[Z, rows] = interp_from_R(R, P, k);
verr = norm(trueprobe - A(:,rows) * (Z * Omega), 'fro') / trueprobenorm;

% Redraw strategy: growing k from the same sketched R does not help (the
% next pivots are close to where R11 becomes ill-conditioned), and an
% unsketched QR would cost O(m n min(m,n)) per block. Instead redraw a dense
% Gaussian sketch (up to two attempts, l then 2l), sharpened by
% escalatepowerits power iterations, and accept it when the probe error is
% below margin*cutoffval.
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

function [R, P] = sketch_qr(A, l, sketcher, q)
% Sketch Y = S*A with l rows (sparse sign by default), optional power
% iterations, then the column-pivoted QR of Y.
m = size(A, 1);
switch lower(sketcher)
    case 'gaussian'
        Y = randn(l, m) * A;
    case 'sparsestack'
        zeta = min(8, min(size(A)));
        Y = hssutil.sparsesign(l, m, zeta) * A;
    case 'sparsestack_slow'
        zeta = min(8, min(size(A)));
        Y = sparsesign_slow(l, m, zeta) * A;
    otherwise
        error("Unknown sketcher type: %s", sketcher);
end
% power iteration (Halko-Martinsson-Tropp Alg. 4.4), re-orthonormalized
% between passes; off by default (q = 0)
for pit = 1:q
    [Qy,~] = qr(Y','econ');
    Z = A * Qy;
    [Qz,~] = qr(Z,'econ');
    Y = Qz' * A;
end
[~, R, P] = qr(Y,'econ','vector');
end

function k = rank_from_R(R, cutoffval)
% numerical rank from the pivoted-QR diagonal, relative to its largest entry
d = abs(diag(R));
if isempty(d) || d(1) == 0
    k = 0;
else
    k = find(d / d(1) >= cutoffval, 1, 'last');
    if isempty(k), k = 0; end
end
end
