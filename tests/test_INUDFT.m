%test_inudft: 
clear all
close all
%%
n  = 2^10+57; 
nd = linspace(0, 2*pi, n+1); 
nd = nd(1:end-1).'; 
nd = nd + (nd(2)-nd(1))*(.5)*rand(n, 1); 
nd = exp(-1i*nd);



%clf
%plot(nd, 'x')
%hold on
%plot(exp(1i*linspace(0, 2*pi, 100)), '.-k');

%%
b = rand(n,1); %rhs 
V = nd.^(0:n-1);
vt = V\b;

v = INUDFT(nd,n,b);
norm(v - vt)
