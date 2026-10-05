function [x, r] = tikhonov_solve(H, b, lam, w, s)
%TIKHONOV_SOLVE  x = argmin ||S (H x - b)||^2 + lam^2 ||W x||^2,  lam > 0,
%   S = diag(s), W = diag(w) (pass [] for identity); H wide, square or tall.
%   r returns the residual S(Hx - b).
%   Thin wrapper of the method @hss/tikhonov,
%   [x, r] = tikhonov(H, b, lam, 'Weight', w, 'DataWeight', s), which also
%   takes leaf blocks for W and S and caches its factors for the last lam.
if nargin < 4, w = []; end
if nargin < 5, s = []; end
[x, r] = tikhonov(H, b, lam, 'Weight', w, 'DataWeight', s);
end
