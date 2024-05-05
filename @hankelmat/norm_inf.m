function h = norm_inf(H)
%NORM_INF Efficiently computes the infinity-norm for Hankel matrices. 
%   The inf-norm is the maximum absolute row sum. 
% 
%   We compute by noting that the infinity norm is invariant under left
%   multiplication by the exchange matrix (resulting in a left-right flip),
%   and so we are able to reuse the inf-norm of Toeplitz!
    h = norm(fliplr(H), 'inf');
end

