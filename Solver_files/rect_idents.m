%rectangular identities
%%
% for the rectangular Toeplitz, we need
% that the diagonalization results in a matrix of the
% form diag(w^2, w^4, ..., w^n), where w = exp(1i*pi/n). 
clear all
close all
%%
n = 5; 
m = 10; 
tc = rand(m, 1); 
tr = rand(n,1); 
tr(1) = tc(1); 
T = toeplitz(tc, tr); 
%%
wn = exp(1i*pi/n);
wm = exp(1i*pi/m);
Z = circshift(eye(n), 1); 
Dr = diag(wn.^(2*(1:n))); %scaling matrix
%%
% (F_n)_jk = 1/sqrt(n) w^2(j-1)(k-1). Apply Fn to v as sqrt(n)*ifft(v). 
% Fn' applies to v as fft(v)/sqrt(n). 

D1 = fft(ifft(Dr*Z*Dr').').';
D1(abs(D1)<1e-8)=0;
spy(D1)
%%
dd = diag(D1);
ddtru = wn.^(2*(1:n));
norm(dd-ddtru.')
%%
% now try the left factor Zg: 
q = exp(1i*2*pi*pi/2/m);
Zg = circshift(eye(m),1); Zg(1,end) = q^(m);
Q = diag(q.^((1:m)));
D2 = fft(ifft(Q*Zg*Q').').';
D2(abs(D2)<1e-8)=0;
spy(D2)
d2 = diag(D2);
trud2 = q*(wm.^(2*(0:m-1))).';
norm(d2-trud2)
%% 
clf
plot(d2, 'kx')
hold on
plot(dd, 'ob')
%%
% now we want to check ZgT-TZ1 = GH'
% get the generators correctly

[G, H] = toep_gens_rect(tr, tc, q^m);
R = Zg*T-T*Z;

norm(R-G*H')




