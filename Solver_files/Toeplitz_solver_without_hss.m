% test rect Toeplitz: 
%%
% test the basic premise of the solver using Cauchy identities
clear all
close all
%% 
% problem setup
n = 2^10; 
m = 2*n; 
tr = randn(n,1); 
tc = randn(m,1); 
tc(1)=tr(1);
Tm = toeplitzmat(tc,tr);
T = full(toeplitz(tc, tr));
xt = randn(n,1);
bt = T*xt; 
%%
%  let's try doing this with our Cauchy structures. 
% the basic idea is that we solve 
% Cy = b, where C is the Cauchy-like matrix (FmQTW'Fn = C)

%% 

wn = exp(pi*1i/n); 
wm = exp(pi*1i/m);
q = exp(1i*2*pi*pi/2/m); % ensures that the nodes are shifted off ROU.

% get generators:
[GG, LL] = toep_gens_rect(tr, tc,q^m); 
        
%transform to Cauchy-like generators:
 Qr = spdiags((wn.^(2*(1:n))).', 0, n,n);
 Ql = spdiags((q.^(1:m)).', 0, m,m);
%%
 G = sqrt(m)*ifft(Ql*GG);
 L = sqrt(n)*ifft(Qr*LL);

 %%
 % construct the Cauchy matrix: 
 nodes = q*(wm.^(2*(0:m-1))).';
 rou = wn.^(2*(1:n)).';

 %%
 % now construct Cauchy matrix: 
 C = (1./(nodes-rou.')).*(G*L'); 
 Dn = diag(nodes);
 Dr = diag(rou);
 RHS = Dn*C-C*Dr; 
 norm(RHS-G*L')
 %%
 Zg = circshift(eye(m),1);
 Zg(1,end)=q^m; 
 testnodes = diag(fft(ifft(Ql*Zg*Ql').').'); 
 norm(testnodes-nodes)
 %%
 % check the matrix transformation: 
 Tm = Ql*T*Qr';
 Ct = (fft(ifft(Tm).').')/sqrt(n)*sqrt(m);
norm(C-Ct)
 %%

 % now modify RHS: 
 b = sqrt(m)*ifft(Ql*bt);
 y = C\b;
 x = Qr'*fft(y)/sqrt(n);
 x = real(x);
 %%

 

