function x = minnorm_dense(A, b)
%MINNORM_DENSE  Backward-stable dense minimum-norm solution (A full row rank):
%   A' = Q R  ->  x = Q (R' \ b).  (Same answer as lsqminnorm / pinv, also in Octave.)
[Q, R] = qr(A', 0);
x = Q * (R' \ b);
end
