function x = mp_weighted_minnorm(H, b, w, Lblocks)
%MP_WEIGHTED_MINNORM  x = argmin ||L x||_2  s.t.  H x = b  (H wide or square, full row rank)
%   L = diag(w)               : mp_weighted_minnorm(H, b, w)
%   L = blkdiag(Lblocks{:})   : mp_weighted_minnorm(H, b, [], Lblocks)
%   Thin wrapper kept for the suite: the method is now @hss/minnorm,
%   x = minnorm(H, b, 'Weight', w) or minnorm(H, b, 'Weight', Lblocks)
%   (memo 2; factors are cached in H for repeated calls with the same weight).
if nargin < 4 || isempty(Lblocks)
    x = minnorm(H, b, 'Weight', w);
else
    x = minnorm(H, b, 'Weight', Lblocks);
end
end
