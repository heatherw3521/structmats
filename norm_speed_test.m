% Levent Batakci 4/22/2024
%
% Just a little speed-test for the 2-norm

n = 1e8;
T = toeplitzmat(rand(n,1));
% A = full(T);

%% Compute norms
h = norm(T)
% h2 = norm(A);