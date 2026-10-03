function [Z, cols] = mp_id_cols(A, tol, kmax)
%MP_ID_COLS  Deterministic column interpolative decomposition
%   A ~= A(:,cols) * Z,  Z (k x n) with Z(:,cols) = I.
%   Rank rule: k = #{j : |R_jj| >= tol*|R_11|} of the column-pivoted QR
%   (the same relative-threshold rule as +hssutil/inter_decompv3.m, but on A
%   itself rather than on a randomized sketch, and without its rank cap).
if nargin < 3, kmax = inf; end
[m, n] = size(A);
if m == 0 || n == 0
    Z = zeros(0, n); cols = zeros(1, 0);
    return
end
[R, P] = mp_pivqr(A);
d = abs(diag(R));
if isempty(d) || d(1) == 0
    k = 0;
else
    k = sum(d / d(1) >= tol);
end
k = min(k, kmax);
if k == 0
    Z = zeros(0, n); cols = zeros(1, 0);
    return
end
T = R(1:k, 1:k) \ R(1:k, k+1:end);
Zp = [eye(k), T];
Z = zeros(k, n, 'like', Zp);
Z(:, P) = Zp;
cols = P(1:k);
end
