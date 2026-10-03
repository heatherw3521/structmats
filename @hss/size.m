function varargout = size(H, varargin)
%SIZE  Size of an HSS matrix, with the calling forms of the built-in size:
%   s = size(H), [m, n] = size(H), m = size(H, 1), s = size(H, [1 2]),
%   [m, n] = size(H, 1, 2). Dimensions beyond 2 are 1.
s = H.sz;
if isempty(s), s = [0 0]; end
if nargin > 1
    dims = [varargin{:}];
    if any(dims < 1 | dims ~= fix(dims))
        error('hss:size:dim', 'Dimension argument must be a positive integer.');
    end
    out = ones(1, numel(dims));
    out(dims <= 2) = s(dims(dims <= 2));
    if nargout <= 1
        varargout{1} = out;
    else
        for k = 1:nargout, varargout{k} = out(k); end
    end
    return
end
if nargout <= 1
    varargout{1} = s;
else
    varargout{1} = s(1);
    varargout{2} = s(2);
    for k = 3:nargout, varargout{k} = 1; end
end
end
