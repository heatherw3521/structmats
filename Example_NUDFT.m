% Examples: Non-uniform DFT matrices 

% Struct solve will eventually be a collection of fast solvers for 
% matrices with special displacement structures. These examples show how to
% use the code as it currently exists.
%
% This example file is for systems involving type II non-uniform 
% discrete Fourier transform matrices.
%%
% Non-uniform type-II DFT matrix: Solve Vx = b, where
% 
%  V is a vandermonde matrix with V(j,k) = exp(-2i*pi*p(j)*(k-1)). 
%
% set up the problem: 
n = 1024; % number of modes
m = 2*n; 
mds = 0:n-1;
p = rand(m,1); %sample locations for non-unif DFT 
V = exp(-2*pi*1i*p.*mds); b = V*rand(n,1); %b = rhs
%%
% solve with structmat: 
x = structsolv_nudft2(p, n, b); 

% check residual: 
norm(V*x-b)./norm(b)
%%
% solve with structmat, add tol parameter: 
x = structsolv_nudft2(p, n, b, 'tol', 1e-6);
norm(V*x-b)/norm(b)

%% Multiple RHS
% To solve VX = B, where B is a collection of multiple right
% hand sides: First, solve for a subset of RHSs) and save the factorization 
% information. 
%
% Here's an example showing how to do that:
XX = rand(n,5);
B = V*XX; 

[L, P, x] = structsolv_nudft2(p, n, B(:,1), 'tol', 1e-6); 
% L is a factorization obj. 
% P is permutation information related to the construction of L. 

%%
% Now that we have the factorization information, we can apply it to 
% additional RHSs with the 'nudft2_solve' command: 

X = nudft2_solve(L,P,B(:,2:end)); 

%%
norm(V*X - B(:,2:end))/norm(B(:,2:end))










