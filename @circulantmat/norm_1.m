function h = norm_1(C)
%TOEP_1_NORM Efficiently computes the 1-norm for Circulant matrices.
%   The 1-norm is the maximum absolute column sum. For circulant matrices,
%   this is the same for every column and row! 
    h = sum(abs(C.tc));
end

