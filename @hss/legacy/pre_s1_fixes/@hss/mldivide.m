function [x, H1] = mldivide(H1, H2)
% backslash: x = H\b for an HSS matrix H and a dense b (one or more columns)
%   square H : the solution of H*x = b
%   wide H   : the minimum-norm solution, argmin ||x|| s.t. H*x = b
% Both use the ULV factorization in @hss/private/hss_ulvminnormsolve.m. The
% factors are kept with H: the first H\b factors, every later H\b with the
% same (unmodified) H only solves, O(n r). clearfactors(H) frees them.
% For weighted minimum norm see minnorm; for regularized/least-squares
% problems (any shape) see tikhonov.

if isa(H1,'hss')
    if isscalar(H2)
        % scalar mult
        error('hss scalar mat multiplication is not yet supported')
        %H = 1/H2* H1;
    elseif isa(H2, 'hss')
        error('hss hss backsolver is not yet supported')
        %H = hss_mldivide(H1, H2);
    else
        if H1.size(1) > H1.size(2)
            error('hss:mldivide:tall', ['overdetermined H\\b is not yet implemented; ' ...
                'for regularized least squares use tikhonov(H, b, lambda).'])
        end
        x = hss_ulvminnormsolve(hss_cachedfactor(H1), H2);
    end
else
    error('hss inverse is not yet supported')
    % H = inv(H2)/inv(H1);
end
end
