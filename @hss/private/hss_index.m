function s = hss_index(H, subs)
%HSS_INDEX  Entries of H for the subscripts of H(...), with MATLAB semantics.
%   H(I,J): I, J may be ':', a logical mask, or positive integers (any order,
%   repeats allowed); further subscripts must be 1. H(K): linear indexing,
%   column-major as for a dense matrix; the result takes the shape of K,
%   or of H's orientation when H is a vector, and H(:) is a column.
sz = H.sz;
if numel(subs) == 1
    k = subs{1};
    N = prod(sz);
    if iscolon(k)
        lin = 1:N; shp = [N 1];
    else
        lin = resolve(k, N);
        if islogical(k)
            shp = [numel(lin) 1];
            if isrow(k), shp = [1 numel(lin)]; end
        else
            shp = size(k);
        end
        if any(sz == 1) && isvector(k)          % vector H: result keeps H's orientation
            if sz(1) == 1, shp = [1 numel(lin)]; else, shp = [numel(lin) 1]; end
        end
    end
    [i, j] = ind2sub(sz, lin);
    [iu, ~, ii] = unique(i); [ju, ~, jj] = unique(j);
    S = block(H, iu, ju);
    s = reshape(S(sub2ind(size(S), ii(:), jj(:))), shp);
    return
end
for d = 3:numel(subs)                            % trailing subscripts: only 1 or ':'
    if ~(iscolon(subs{d}) || isequal(subs{d}, 1) || isequal(subs{d}, true))
        error('hss:index:outOfRange', 'Index in position %d exceeds array bounds (must not exceed 1).', d);
    end
end
rows = resolve(subs{1}, sz(1));
cols = resolve(subs{2}, sz(2));
s = block(H, rows, cols);
end

function S = block(H, rows, cols)
if isempty(rows) || isempty(cols)
    S = zeros(numel(rows), numel(cols));
else
    S = extract(H, rows, cols, 0, 0);
end
end

function tf = iscolon(x)
tf = (ischar(x) || isstring(x)) && isequal(char(x), ':');
end

function idx = resolve(x, n)
% subscript -> row vector of integers in 1..n, or an error
if iscolon(x)
    idx = 1:n;
elseif islogical(x)
    if numel(x) > n && any(x(n+1:end))
        error('hss:index:outOfRange', 'The logical index contains a true value outside the array bounds (%d).', n);
    end
    idx = find(x(:)).';
elseif isnumeric(x)
    idx = double(x(:)).';
    if any(~isreal(x)) || any(idx ~= fix(idx)) || any(~isfinite(idx)) || any(idx < 1)
        error('hss:index:invalid', 'Array indices must be positive integers or logical values.');
    end
    if any(idx > n)
        error('hss:index:outOfRange', 'Index exceeds the number of rows/columns (%d).', n);
    end
else
    error('hss:index:invalid', 'Array indices must be positive integers or logical values.');
end
end
