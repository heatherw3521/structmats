function h = norm_fro(C)
%NORM_FRO Efficiently computes the frobenius-norm for Circulant matrices.
%   The Frobenius norm is defined as the square-root of the sum of 
%   the squared entries.
    h = sqrt(size(C,1) * sum((C.tc).^ 2)); % Each element appears m times!
end

