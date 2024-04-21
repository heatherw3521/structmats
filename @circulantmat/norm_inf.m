function h = norm_inf(C)
%NORM_INF Efficiently computes the infinity-norm for Circulant matrices.
%   The 1-norm is the maximum absolute row sum. For circulant matrices,
%   this is the same for every row and column! 
    h = sum(abs(C.tc));
end

