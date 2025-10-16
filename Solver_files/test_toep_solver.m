% test fast Toeplitz solver: 
clear all
close all
%%
n = 2^12; 
m = 2^12; 
%m = n;
tr = randn(n,1) + 1i*randn(n,1); 
tc = randn(m,1); 
tc(1)=tr(1);
Tm = toeplitzmat(tc,tr);
T = full(toeplitz(tc, tr));
xt = randn(n,1)+1i*randn(n,1);
b = T*xt; 
%%
tt = tic;
x = structsolv_toeplitz(tc,tr, b);
time_solve = toc(tt)
%x = Tm\b;
%%
norm(xt-x)/norm(xt)
%%