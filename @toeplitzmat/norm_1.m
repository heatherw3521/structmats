function h = norm_1(T)
%NORM_1 Efficiently computes the 1-norm for Toeplitz matrices.
%   The 1-norm is the maximum absolute column sum. 
%
%   We compute this by noting that the the maximum absolute column sum is
%   the maximum absolute row sum (infinity-norm!) of the transpose.

    h = norm(T.', 'inf');
end

