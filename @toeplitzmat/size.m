function s = size(T, dim)
%SIZE for TOEPLITZMAT objects
%  returns the size of toeplitz matrix

    s = [size(T.tc, 1), size(T.tr, 2)];
    if(exist('dim','var'))
        if(dim <= 0 || dim >= 3)
            error(['TOEPLITZMAT:size:baddimension', ...
                   'Cannot get size of toeplitzmatrix in dimension %s ',""+dim]);
        end
       s = s(dim);
    end
    
end