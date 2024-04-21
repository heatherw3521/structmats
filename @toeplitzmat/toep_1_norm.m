function h = toep_1_norm(T)
%TOEP_1_NORM Quickly computes the 1-norm for Toeplitz matrices.
%   The 1-norm is the maximum absolute column sum. 
% 
%   For toeplitz matrices, the 1-norm and infinity-norm are the same! This
%   can most easily be seen by the fact that ||A||_1 = ||A.'||_inf, and the
%   symmetry in the inf-norm algorithm, which is invariant under swapping
%   T.tc and T.tr.

    h = toep_inf_norm(T);
end

