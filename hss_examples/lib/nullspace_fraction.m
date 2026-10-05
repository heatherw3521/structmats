function f = nullspace_fraction(A, x)
%NULLSPACE_FRACTION  ||P_null(A) x|| / ||x|| (dense): 0 for an exact min-norm solution.
[Q, ~] = qr(A', 0);
f = norm(x - Q*(Q'*x)) / norm(x);
end
