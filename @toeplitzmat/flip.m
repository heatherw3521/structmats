function h = flip(T)
%FLIP up-down flip for TOEPLITZMAT
%   The flip of a Toeplitz matrix is a Hankel matrix!
    
h = hankelmat(flip(T.tc), T.tr);
end

