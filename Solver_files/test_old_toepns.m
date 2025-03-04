%test toep sizes
clear all
close all
tol = 1e-10; 
m = 1028; 
tc = rand(m,1);
tr = rand(m,1);
tc(1) = tr(1);
b = rand(m,1); 
xt = full(toeplitz(tc,tr))\b;
[xs, ~,solvt] = Toeplitz_solve_ns(tc,tr,b, 'tol', tol);
norm(xs-xt)
%%