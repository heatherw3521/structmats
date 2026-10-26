function S = sparsesign(d, m, zeta)
%SPARSESIGN  Sparse sign embedding: a d x m sparse matrix whose every column
%   has zeta nonzeros, +-1/sqrt(zeta) with random signs, in distinct random
%   rows (zeta is capped at d).
%
%   Portable MATLAB version of the compiled sparsesign (sparsesign.c). Where
%   a compiled version exists for the platform (e.g. sparsesign.mexmaca64),
%   MATLAB uses it instead of this file; elsewhere this file is used, with
%   the same distribution. To compile for your platform: mex sparsesign.c
zeta = min(zeta, d);
if zeta <= 0 || m == 0
    S = sparse(d, m);
    return
end
if d <= 4*zeta
    % few rows: a random permutation of the rows of every column
    [~, idx] = sort(rand(d, m), 1);
    rows = idx(1:zeta, :);
else
    % many rows: draw with replacement and redraw the columns with repeats
    rows = randi(d, zeta, m);
    while true
        dup = any(diff(sort(rows, 1), 1, 1) == 0, 1);
        if ~any(dup), break; end
        rows(:, dup) = randi(d, zeta, nnz(dup));
    end
end
cols = repmat(1:m, zeta, 1);
vals = (2*randi(2, zeta, m) - 3) / sqrt(zeta);
S = sparse(rows(:), cols(:), vals(:), d, m);
end
