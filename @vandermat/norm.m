function h = norm(V,varargin)
%NORM Norm function for Vandermonde matrices. For certain special cases, e.g.
%    the infinity-norm and Frobenius norm, we can compute the
%    value very efficiently.

    if(isempty(varargin))
        h = builtin('norm',full(V));
    else
        switch lower(varargin{1})
            case {'fro', 'frobenius'}
                h = norm_fro(V);
            case {'inf'}
                h = norm_inf(V);
            otherwise
                h = builtin('norm',full(V),varargin{:});
        end
    end

end

