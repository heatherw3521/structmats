% Example_Toeplitz: 
% A brief example showing how to use the fast Toeplitz solver:
%
% set up a Toeplitz matrix:
clear all
%n = 2^12;   
n = 1e4; 
tr = rand(n,1) + 1i*rand(n,1);
tc = rand(n,1) + 1i*rand(n,1); tc(1)= tr(1); 
T = toeplitz(tc, tr); 
xt = rand(n,1);
b = T*xt; 
 
%%
% try the solver: 
xns = structsolv_toeplitz(tc, tr,b);
norm(xt-xns)/norm(xt)
%% 
% You can also set a desired tolerance (The default is 1e-10). 
xns = structsolv_toeplitz(tc, tr,b,'tol',1e-5);
norm(xt-xns)/norm(xt)


%% Multiple RHS
% For multiple RHS, we can save the factorization information after one solve, 
% and then apply the factorization directly to additional RHSs. First, 
% solve with one RHS:
Xt = rand(n, 5); 
B = T*Xt; 
[L,xns] = structsolv_toeplitz(tc, tr, B(:,1), 'tol', 1e-5); 
%%
% Now, we can apply to other RHSs with the 'toeplitz_solve' command: 
X = toeplitz_solve(L, B); 
norm(Xt-X)/norm(Xt)




