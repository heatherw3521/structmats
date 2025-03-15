% test rect Toeplitz: 
clear all
close all
%%
n = 2^10; 
m = 2*n; 
tr = randn(n,1); 
tc = randn(m,1); 
tc(1)=tr(1);
Tm = toeplitzmat(tc,tr);
T = full(toeplitz(tc, tr));
xt = randn(n,1);
b = T*xt; 
%%
x = structsolv_toeplitz(tc,tr, b);
%x = Tm\b;
%%
norm(xt-x)/norm(xt)
%%