function [nrm, its] = normest(H, tol, maxit)
%NORMEST  Estimate of the 2-norm of an HSS matrix by power iteration on
%   H'*H, using only matvecs (O(n r) each).
%   nrm = normest(H) stops when two successive estimates agree to 1e-6;
%   normest(H, tol, maxit) sets the relative tolerance and the maximum
%   number of iterations (default 100). [nrm, its] also returns the number
%   of iterations.
if nargin < 2 || isempty(tol), tol = 1e-6; end
if nargin < 3 || isempty(maxit), maxit = 100; end
hss_assertroot(H, 'normest');
n = H.sz(2);
Ht = ctranspose(H);
x = randn(n, 1);
x = x / norm(x);
nrm = 0;
for its = 1:maxit
    y = hss_matvec(H, x);
    ny = norm(y);
    if ny == 0
        nrm = 0;
        return
    end
    x = hss_matvec(Ht, y);
    nrm0 = nrm;
    nrm = norm(x) / ny;              % Rayleigh-quotient estimate of ||H||
    x = x / norm(x);
    if abs(nrm - nrm0) <= tol * nrm
        return
    end
end
end
