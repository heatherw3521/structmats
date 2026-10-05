function [x, H1] = mldivide(H1, H2)
%MLDIVIDE  x = H\b for an HSS matrix H and a dense b (one or more columns).
%   square H : the solution of H*x = b
%   wide H   : the minimum-norm solution, argmin ||x|| s.t. H*x = b
%   Both use a ULV factorization (private/hss_ulvminnormsolve.m). The factors
%   are kept with H: the first H\b factors, every later H\b with the same
%   (unmodified) H only solves, in O(n r). clearfactors(H) frees them.
%   For weighted minimum norm see minnorm; for regularized least squares
%   (any shape) see tikhonov.

if isa(H1,'hss')
    hss_assertroot(H1, 'H\b');
    if isscalar(H2)
        error('hss:mldivide:scalar', 'H\\s with a scalar s is not supported.')
    elseif isa(H2, 'hss')
        error('hss:mldivide:hss', 'H1\\H2 with two HSS matrices is not supported.')
    else
        if H1.sz(1) > H1.sz(2)
            error('hss:mldivide:tall', ['overdetermined H\\b is not yet implemented; ' ...
                'for regularized least squares use tikhonov(H, b, lambda).'])
        end
        x = hss_ulvminnormsolve(hss_cachedfactor(H1), H2);
    end
else
    error('hss:mldivide:dense', 'A\\H with a non-HSS A is not supported.')
end
end
