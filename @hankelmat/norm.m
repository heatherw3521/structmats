function h = norm(H,varargin)
%NORM Norm function for Hankel matrices. For certain special cases, e.g.
%    the 1-norm, infinity-norm, and Frobenius norm, we can compute the
%    value very efficiently.

    if(isempty(varargin))
        h = builtin('norm',full(H));
    else
        switch lower(varargin{1})
            case {'fro', 'frobenius'}
                h = hank_frobenius_norm(H);
            case {'1'}
                h = hank_1_norm(H);
            case {'inf'}
                h = hank_inf_norm(H);
            otherwise
                h = builtin('norm',full(H),varargin{:});
        end
    end

end

