clear all
close all
%%
%  Zm*T-T*Zn = GH^*
n = 10; 
m = 20; 
t1 = rand(m, 1); 
t2 = rand(n,1); 
t2(1) = t1(1); 
T = toeplitz(t1, t2); 
Zm = circshift(eye(m), 1); Zm(1,end) = -1;
Zn = circshift(eye(n),1);  
%%
% 
% check transforms on Zn and Zm: 
wn = exp(1i*pi/n); 
wm = exp(1i*pi/m);
D1t = wn.^(2*(0:n-1)); 
Dmt = wm.^(2*(1:m)-1); 
D0 = diag(wm.^(0:m-1)); 

%%
%
% First Zm
Dm = D0*Zm*D0'; 
Dm = fft((ifft(Dm)).').'; 
norm(diag(Dm) - Dmt.')
%%
%
% Now Zn
D1 = fft((ifft(Zn)).').'; 
norm(diag(D1) - D1t.')
%%
% 