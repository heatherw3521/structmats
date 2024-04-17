function H = nearby_hankel(A)
%NEARBY_HANKEL Finds the closest Hankel matrix to A, in the sense of
%   the Frobenius norm. It can be shown that the closest matrix is obtained
%   by simply averaging along the diagonals!
%
%   TODO: Add proof to devlog. It's just basic calculus.
%
%   We do this by reusing nearby_toeplitz and the fact that the flip of a
%   Hankel matrix is Toeplitz.

H = flip(nearby_toeplitz(flip(A)));
end

