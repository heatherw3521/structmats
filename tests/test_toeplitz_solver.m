% basic Toeplitz solver tests

clear all
n = 2^10 ;   % pick non power of 2
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

%%
% try the nonsingular version: 
xns = structsolv_toeplitz(tc, tr,b,'tol',tol );
%xns = Toeplitz_solve_ns(tc,tr,b, 'tol', tol);

norm(x-xns)/norm(x)



