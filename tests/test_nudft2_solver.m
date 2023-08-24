%test_inudft: 
clear all
close all
%%
n  = 2^10+57; 
nd = linspace(0, 1, n+1); 
nd = nd(1:end-1).'; 
nd = nd + (nd(2)-nd(1))*(.5)*rand(n, 1); 
ndv = exp(-2i*pi*nd);

%clf
%plot(ndv, 'x')
%hold on
%plot(exp(1i*linspace(0, 2*pi, 100)), '.-k');

%%
b = rand(n,1); %rhs 
V = ndv.^(0:n-1);
vt = V\b;

v = structsolv_nudft2(nd,n,b);
norm(v - vt)
