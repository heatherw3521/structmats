function h = norm(T,varargin)
%NORM Norm function for Toeplitz matrices. For certain special cases, e.g.
%    the 1-norm, infinity-norm, and Frobenius norm, we can compute the
%    value very efficiently.

    if(isempty(varargin))
        h = builtin('norm',full(T));
    else
        switch lower(varargin{1})
            case {'fro', 'frobenius'}
                h = toep_frobenius_norm(T);
            case {'1'}
                h = toep_1_norm(T);
            case {'inf'}
                h = toep_inf_norm(T);
            otherwise
                h = builtin('norm',full(T),varargin);
        end
    end

end

