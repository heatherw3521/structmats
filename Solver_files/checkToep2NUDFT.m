% check the identities involved in converting the Toeplitz Cauchy-like
% matrix to one that can be passed into the NUDFT solver. 

%%
clear all 
close all
%%
n = 3; 
m = 6; % when this is double n then certain nodes and rou coincide. 
% this messes up the Sylvester equation identity. 
tr = rand(n, 1); 
tc = rand(m,1); 
tc(1) = tr(1); 
T = toeplitz(tc, tr); 
x = rand(n,1);
b = T*x; 

%%
% check disp struct:

%  Zm*T-T*Zn = GH^*
Zn = circshift(eye(n), 1); Zn(1,end) = -1;
Zm = circshift(eye(m),1);


%%
% now modify to Cauchy generators: 
wn = exp(1i*pi/n); 
wm = exp(1i*pi/m);
D1 = wm.^(2*(0:m-1));
Dm1 = wn.^(2*(1:n)-1);
D0 = diag(wn.^(0:n-1));
Dnodes = wn*diag(D1);
Dr = wn*diag(Dm1);
%%
[G, H] = toep_gens_rect(tr, tc); 
cG = wn*ifft(G)*sqrt(m);
cH = ifft(D0*H)*sqrt(n);
%%
% check that this works: 
C = (fft((sqrt(m)*ifft(T)*D0').').')./sqrt(n);

Cdisp = Dnodes*C-C*Dr;
%%
norm(Cdisp-cG*cH')
%%
Ct = (1./(diag(Dnodes)-diag(Dr).')).*(cG*cH'); 
%%
norm(Ct-C)
%%
bhat = sqrt(m)*ifft(b);
xhat = C\bhat; 
%%
xc = D0'*fft(xhat)/sqrt(n);
norm(real(xc)-x)
%%
% now check if it works with permutation: 
 %permute so that nodes are ordered wrt argument:
 nodes = diag(Dnodes);
 args = mod(angle(nodes), 2*pi);
 [args, p] = sort(args); 
%circshift so nodes are ordered appropriately for subdivisions. 
 kk = find(args > pi/n, 1); 
 p = circshift(p, -kk+1); 
 pnodes = nodes(p); % D = diag(nodes)
 cGp = cG(p,:);
 bhatm = bhat(p,:);
 Dnodesp = diag(pnodes);

 %%
 % check disp: 
 Cp  = C(p,:);
 pdisp = Dnodesp*Cp-Cp*Dr;
 norm(pdisp-cGp*cH')
 %%
bhatp = bhat(p,:);
hatxp = Cp\bhatp; 
xp = D0'*fft(hatxp)/sqrt(n);



