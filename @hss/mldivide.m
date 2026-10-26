function [x, H1] = mldivide(H1, H2)
%MLDIVIDE  x = H\b for an HSS matrix H and a dense b (one or more columns).
%   square H : the solution of H*x = b
%   wide H   : the minimum-norm solution, argmin ||x|| s.t. H*x = b
%   tall H   : the least-squares solution, argmin ||H*x - b|| (H of full
%              column rank), as the adjoint of the minimum-norm solve with H'
%   All use a ULV factorization (private/hss_ulvminnormsolve.m). The factors
%   are kept with H: the first H\b factors, every later H\b with the same
%   (unmodified) H only solves, in O(n r). clearfactors(H) frees them.
%   For weighted minimum norm see minnorm; for regularized least squares
%   (any shape) see tikhonov.

if isa(H1,'hss')
    hss_assertroot(H1, 'H\b');
    if isscalar(H2) && H1.sz(1) ~= 1
        error('hss:mldivide:scalar', 'H\\s with a scalar s is not supported.')
    elseif isa(H2, 'hss')
        error('hss:mldivide:hss', 'H1\\H2 with two HSS matrices is not supported.')
    else
        if size(H2, 1) ~= H1.sz(1)
            error('hss:mldivide:dimension', 'b has %d rows; H has %d rows.', size(H2, 1), H1.sz(1));
        end
        if H1.sz(1) > H1.sz(2)
            % least squares: x = pinv(H)*b = pinv(H')'*b
            msg = ['H does not have full column rank, so the least-squares solution ' ...
                'is not unique; use tikhonov(H, b, lambda).'];
            try
                F = hss_cachedadjoint(H1);
            catch err
                if strcmp(err.identifier, 'hss_ulvminnormsolve:rankDeficient')
                    error('hss:mldivide:rankDeficient', msg)
                end
                rethrow(err)
            end
            if ~hss_ulv_fullrank(F)
                error('hss:mldivide:rankDeficient', msg)
            end
            x = hss_ulv_adjoint(F, H2);
        else
            x = hss_ulvminnormsolve(hss_cachedfactor(H1), H2);
        end
    end
else
    error('hss:mldivide:dense', 'A\\H with a non-HSS A is not supported.')
end
end
