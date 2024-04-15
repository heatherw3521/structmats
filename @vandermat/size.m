function s = size(V, dim)
%SIZE for VANDERMAT objects
%  returns the size of vandermonde matrix

    s = [size(V.x, 1), V.n];
    if(exist('dim','var'))
        if(dim <= 0 || dim >= 3)
            error(['VANDERMAT:size:baddimension', ...
                   'Cannot get size of vandermat in dimension %s '],""+dim);
        end
       s = s(dim);
    end
    
end