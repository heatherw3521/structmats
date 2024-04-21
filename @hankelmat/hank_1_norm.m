function h = hank_1_norm(H)
%HANK_1_NORM Quickly computes the 1-norm for Hankel matrices.
%   The 1-norm is the maximum absolute column sum. 
%
%   We compute by noting that the 1 norm is invariant under right
%   multiplication by the exchange matrix (resulting in a up-down flip),
%   and so we are able to reuse the inf-norm of Toeplitz!
    h = norm(flip(H), 1);
end

