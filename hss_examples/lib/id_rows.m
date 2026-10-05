function [Z, rows] = id_rows(A, tol, kmax)
%ID_ROWS  Row interpolative decomposition A ~= Z * A(rows,:), Z(rows,:) = I.
if nargin < 3, kmax = inf; end
[Zt, rows] = id_cols(A', tol, kmax);
Z = Zt';
end
