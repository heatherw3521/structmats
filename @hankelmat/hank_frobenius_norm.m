function h = hank_frobenius_norm(H)
%HANK_FROBENIUS_NORM Efficiently computes the frobenius norm of a Hankel
%   matrix. The Frobenius norm is defined as the square-root of the sum of 
%   the squared entries.
%
%   We compute by noting that the frobenius norm is invariant under
%   permutation of the elements. Thus, we can multiply by the exchange 
%   matrix (resulting in a flip), and reuse the frobenius-norm of Toeplitz!
    h = norm(flip(H), 'fro');
end

