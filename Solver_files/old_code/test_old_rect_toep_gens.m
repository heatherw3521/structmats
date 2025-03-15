clear all
close all
%%
m = 48; 
n = 40; 
Z1 = [zeros(1,m-1); eye(m-1)]; 
Z1 = [Z1, [1; zeros(m-1,1)]];
Z2 = [zeros(1,n-1); eye(n-1)]; 
Z2 = [Z2, [-1; zeros(n-1,1)]];
%%
tc = rand(m,1); 
tr = rand(n, 1); 
tc(1)= tr(1);
T = full(toeplitz(tc,tr));
TrueRHS = Z1*T-T*Z2; 
[G, H] = toep_gens_rect(tr, tc);
norm(G*H'-TrueRHS)
%%
