% basic Toeplitz solver test
% for now, the Toeplitz matrix must be square
clear all
n = 2^10 +47 ;   
testcomplex = 1;     
tr = rand(n,1);
tc = rand(n,1);
if testcomplex, tr = tr + 1i*rand(n,1); tc = tc + 1i*rand(n,1); end
tc(1) = tr(1); 
T = toeplitz(tc, tr); 
x = rand(n,1);
if testcomplex, x = x + 1i*rand(n,1); end
b = T*x; 
tol = 1e-10; 

%%
% try the solver: 
xns = structsolv_toeplitz(tc, tr,b,'tol',tol );
norm(x-xns)/norm(x)

%% try multiple RHS: 

Xt = rand(n, 5); 
B = T*Xt; 
%%
[L,~] = structsolv_toeplitz(tc, tr, B(:,1), 'tol', tol); 
X = Toep_solve(L, B); 
norm(Xt-X)




