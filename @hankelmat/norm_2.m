function h = norm_2(H)
%NORM_2 Efficiently computes the 2-norm for Hankel matrices.
%   The 2-norm is equal to the maximum singular value.
%
%   We do this by using the fact that the 2-norm is invariant under
%   orthogonal transforms, and flipping is one such transform! So, we just
%   use the algorithm from Toeplitz.

    h = norm(flip(H));
end

