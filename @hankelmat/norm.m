function h = norm(H,varargin)
%NORM Norm function for Hankel matrices. For certain special cases, e.g.
%    the 1-norm, infinity-norm, and Frobenius norm, we can compute the
%    value very efficiently.

    if(isempty(varargin))
        h = norm_2(H);
    else
        switch lower(varargin{1})
            case {'fro', 'frobenius'}
                h = norm_fro(H);
            case {'1',1}
                h = norm_1(H);
            case {'inf'}
                h = norm_inf(H);
            case {'2',2}
                h = norm_2(H);
            otherwise
                h = builtin('norm',full(H),varargin{:});
        end
    end

end

