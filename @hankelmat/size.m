function s = size(H, dim)
%SIZE for HANKELMAT objects
%  returns the size of Hankel matrix

    s = [size(H.hc, 1), size(H.hr, 2)];
    if(exist('dim','var'))
        if(dim <= 0 || dim >= 3)
            error( 'TOEPLITZMAT:size:baddimension', ...
                ['Cannot get size of toeplitzmatrix in dimension %s '],""+dim);
        end
       s = s(dim);
    end
    
end