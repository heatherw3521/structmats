function S = nearby(matrixtype,A)
%NEARBY gets the closest Toeplitz/Hankel/etc. matrix to a given input A, in
%   the sense of the Frobenius norm.
%   
%   matrixtype: This is the type of special matrix we will look for
%   A: This is the matrix for which we find a nearby special matrix


switch lower(matrixtype)
    case {'toeplitz','toeplitzmat','t'}
        S = nearby_toeplitz(A);
    case {'hankel','hankelmat','h'}
        S = nearby_hankel(A);
    otherwise
        error(['Currently, the NEARBY function does not support ',...
            'projection to matrices of type %s'], matrixtype)
end
end

