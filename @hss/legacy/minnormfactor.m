function [solve, F] = minnormfactor(H)
%MINNORMFACTOR  Factor a wide (or square) HSS matrix once for minimum-norm solves.
%   solve = minnormfactor(H) factors H with the recursive ULV minimum-norm
%   solver (@hss/private/hss_ulvminnormsolve.m) and returns a function handle:
%   x = solve(b) is the minimum-norm solution of H*x = b, i.e. the same as
%   H\b, but every call after the factorization costs only O(n r) (triangular
%   solves, one HSS product per level, small unitary products). b may have
%   several columns. Use it for repeated solves with one matrix: Douglas-
%   Rachford or ADMM projections, multiple right-hand sides.
%
%   [solve, F] = minnormfactor(H) also returns the stored factors (a struct,
%   one entry per recursion level in F, F.next, F.next.next, ...).
%
%   Example:
%       solve = minnormfactor(H);
%       for k = 1:K, v = v + solve(b - H*v); end    % repeated projections
if H.size(1) > H.size(2)
    error('hss:minnormfactor:tall', ...
        'minnormfactor needs a wide or square H (got %d x %d).', H.size(1), H.size(2));
end
F = hss_ulvminnormsolve(H);
solver = @hss_ulvminnormsolve;      % handle to the private function, bound here
solve = @(b) solver(F, b);
end
