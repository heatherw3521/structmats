function [Z, rows, varargout] = id_librla(A, varargin)
%ID_LIBRLA  Interpolative decomposition (ID) computed with librla.
%
% A ≈ Z*A(rows,:)                 if orientation = 'rows'
% A ≈ A(:,rows)*Z                 if orientation = 'columns'
% A ≈ Z*A(rows,cols)*Z2           if orientation = 'double-sided'
%
% Same interface as hssutil.inter_decompv3, with the sketching done by
% librla (https://github.com/agillman20/librla, which must be on the path).
%
% NAME-VALUE INPUTS:
%   ctype       : 'threshold' (relative tolerance) or 'k' (fixed rank);
%                 default 'threshold'
%   cval        : the tolerance or the rank (default: 1e-12)
%   orientation : 'rows', 'columns', or 'double-sided'
%   variant     : 'id' (default) or 'range'
%       'id'    librla.qr_sketch: range finder, projection of A onto the
%               range basis, pivoted QR of the projection (as in
%               librla.id_sketch).
%       'range' librla.orth_sketch only (one pass over A): the ID is taken
%               from the leading part of the range basis, then checked
%               against A on a random probe vector; if the probe error
%               exceeds escalatemargin*cval, the block is redone with
%               'id' and escalatepowerits power iterations.
%   powerits    : power iterations for librla (default 0)
%   escalatemargin, escalatepowerits : see 'range' (defaults 10 and 1)
%
% In 'threshold' mode the sketch size is chosen here rather than by
% librla's own adaptive loop: librla's loop grows the sketch 4x per step
% and, once that passes the short side of the block, switches to a full
% pivoted QR of A. HSS block rows are short and very wide, often with a
% rank close to the short side, so that fallback is slow. Here the target
% rank starts at 30 and doubles, and at the short side the sketch simply
% has that many columns.
%
% OUTPUTS:
%   Z, rows
%   varargout : {cols, Z2} for the double-sided case

opts = struct('ctype', 'threshold', 'cval', 1e-12, 'orientation', 0, ...
    'variant', 'id', 'powerits', 0, 'escalatemargin', 10, 'escalatepowerits', 1);
if mod(numel(varargin), 2) ~= 0
    error('hssutil:id_librla:args', 'Options must be name-value pairs.');
end
for i = 1:2:numel(varargin)
    name = char(varargin{i});
    if ~isfield(opts, name)
        error('hssutil:id_librla:args', 'Unknown option ''%s''.', name);
    end
    opts.(name) = varargin{i+1};
end
if ~(strcmp(opts.ctype, 'threshold') || strcmp(opts.ctype, 'k'))
    error('hssutil:id_librla:ctype', 'ctype must be ''threshold'' or ''k'' (got ''%s'').', char(opts.ctype));
end
if ~(strcmp(opts.variant, 'id') || strcmp(opts.variant, 'range'))
    error('hssutil:id_librla:variant', 'variant must be ''id'' or ''range'' (got ''%s'').', char(opts.variant));
end

side = opts.orientation;
if isequal(side, 0)
    if size(A,2) >= size(A,1), side = 'rows'; else, side = 'columns'; end
end

switch side
    case 'columns'
        [Z, rows] = column_id(A, opts);
    case 'rows'
        % a row ID of A is a column ID of A.'
        [Zt, rows] = column_id(A.', opts);
        Z = Zt.';
    case 'double-sided'
        [Z, rows] = hssutil.id_librla(A, varargin{:}, 'orientation', 'rows');
        [Z2, cols] = hssutil.id_librla(A(rows,:), varargin{:}, 'orientation', 'columns');
        varargout{1} = cols;
        varargout{2} = Z2;
    otherwise
        error('hssutil:id_librla:orientation', "orientation must be 'rows', 'columns', or 'double-sided'.")
end
end


function [Z, cols] = column_id(A, opts)
% A ≈ A(:,cols)*Z
[m, n] = size(A);
if m == 0 || n == 0 || ~any(A(:)) || (strcmp(opts.ctype, 'k') && opts.cval == 0)
    Z = zeros(0, n, 'like', A); cols = [];
    return
end

if strcmp(opts.ctype, 'k')
    [Z, cols] = qr_id(A, 0, min(opts.cval, min(m, n)), opts.powerits);
elseif strcmp(opts.variant, 'range')
    [Z, cols, ok] = range_id(A, opts);
    if ~ok
        [Z, cols] = qr_id(A, opts.cval, 0, opts.escalatepowerits);
    end
else
    [Z, cols] = qr_id(A, opts.cval, 0, opts.powerits);
end

if numel(cols) == min(m, n)
    warning('hssutil:notLowRank', 'Off-diagonal block may not be low-rank (numerical rank = min(m,n)).')
end
end


function [Z, cols] = qr_id(A, tol, krank, powerits)
% Column ID from librla.qr_sketch (projected pivoted QR) in its rank mode,
% A(:,p) ≈ Q*R. With krank > 0 the rank is fixed; otherwise (tol > 0) the
% target rank doubles until the tolerance is met within it.
s = min(size(A));
if krank > 0
    [R, p] = projected_qr(A, krank, powerits);
    k = krank;
else
    kmax = min(30, s);
    while true
        [R, p] = projected_qr(A, kmax, powerits);
        d = abs(diag(R));
        if isempty(d) || d(1) == 0
            k = 0;
        else
            k = sum(d >= tol * d(1));
        end
        if k < kmax || kmax == s
            break
        end
        kmax = min(2*kmax, s);
    end
end
k = min(k, size(R, 1));
cols = p(1:k);
Z = zeros(k, size(A, 2), 'like', R);
Z(:, p(1:k)) = eye(k);
if k > 0 && k < size(R, 2)
    R11 = R(1:k, 1:k);
    d = abs(diag(R11));
    if min(d) > eps(class(R)) * max(d)
        Z(:, p(k+1:end)) = R11 \ R(1:k, k+1:end);
    else
        Z(:, p(k+1:end)) = pinv(R11) * R(1:k, k+1:end);
    end
end
end


function [R, p] = projected_qr(A, k, powerits)
% librla.qr_sketch in rank mode, with the oversampling capped so that the
% sketch never has more columns than the short side of A
s = min(size(A));
extra = max(0, min(12, s - k));
[~, R, p] = librla.qr_sketch(A, k, 'power_iter', powerits, 'extra_samples', extra);
end


function [Z, cols, ok] = range_id(A, opts)
% ID from librla's basis of the row space of A (columns of A.'): Q has
% orthonormal columns with A ≈ A*Q*Q'. Pick k columns J by a pivoted QR of
% Q(:,1:k)' and set Z = (Q(J,1:k)')^{-1} Q(:,1:k)', so Z(:,J) = I.
[~, n] = size(A);
s = min(size(A));
tol = opts.cval;
kmax = min(30, s);
while true
    extra = max(0, min(12, s - kmax));
    [Q, ~, diagR] = librla.orth_sketch(A.', kmax, 'power_iter', opts.powerits, 'extra_samples', extra);
    if isempty(diagR) || diagR(1) == 0
        k = 0;
    else
        k = sum(diagR >= tol * diagR(1));
    end
    if k <= kmax || kmax == s
        break
    end
    kmax = min(2*kmax, s);
end
k = min([k, kmax, size(Q, 2)]);
ok = false;
if k == 0
    Z = zeros(0, n, 'like', A); cols = []; ok = true;
    return
end
Qk = conj(Q(:, 1:k));      % basis of the row space of A (range of A')
[~, ~, p] = qr(Qk', 0);
cols = p(1:k);
Z = Qk(cols, :)' \ Qk';
Z(:, cols) = eye(k);

% probe check against A itself
w = randn(n, 1);
if ~isreal(A), w = w + 1i*randn(n, 1); end
ref = A * w;
nref = norm(ref);
if nref == 0
    ok = true;
    return
end
err = norm(ref - A(:, cols) * (Z * w)) / nref;
ok = err <= opts.escalatemargin * tol;
end
