function x = weighted_minnorm_solve(H, b, w, Lblocks)
%WEIGHTED_MINNORM_SOLVE  x = argmin ||L x||_2  s.t.  H x = b  (H wide or square, full row rank)
%   L = diag(w)               : weighted_minnorm_solve(H, b, w)
%   L = blkdiag(Lblocks{:})   : weighted_minnorm_solve(H, b, [], Lblocks)
%   Thin wrapper of the method @hss/minnorm,
%   x = minnorm(H, b, 'Weight', w) or minnorm(H, b, 'Weight', Lblocks)
%   (factors are cached in H for repeated calls with the same weight).
if nargin < 4 || isempty(Lblocks)
    x = minnorm(H, b, 'Weight', w);
else
    x = minnorm(H, b, 'Weight', Lblocks);
end
end
