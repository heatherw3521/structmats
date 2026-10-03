function [Z, rows] = mp_id_rows(A, tol, kmax)
%MP_ID_ROWS  Row interpolative decomposition A ~= Z * A(rows,:), Z(rows,:) = I.
if nargin < 3, kmax = inf; end
[Zt, rows] = mp_id_cols(A', tol, kmax);
Z = Zt';
end
