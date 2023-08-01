cd %test multiple RHS: 

clear all
close all
%%
n  = 2^10; m = 2*n;
nd = exp(-1i*2*pi*rand(1,m)).';
%%
N = 20; 
B = rand(m,N); %rhs 
V = nd.^(0:n-1);
Xt = V\B;
s = tic;
[L,p,x] = INUDFT(nd,n,B(:,1));
%%
X = INUDFT_solve(L, p, B);
toc(s)
norm(X- Xt)./norm(Xt)