%square test: this is  a simple test to see if the constructor does 
% the right thing for a two-level tree. 

clear all
close all
%%
m = 32; 
n = 32; 
b = 8; % blocksize
r = 4; %rank
A = rand(m,r)*rand(r,n) + eye(m);
% top half matrix:
% A(1:2*b, :) = rand(2*b, r)*rand(r, n); 
% A(1:b, 1:b) = rand(b,b);
% A(b+1:2*b, b+1:2*b) = rand(b,b); 
% % bottom half matrix: 
% A(2*b+1:end, :) = rand(2*b, r)*rand(r, n); 
% A(2*b+1:3*b, 2*b+1:3*b) = rand(b,b);
% A(3*b+1:4*b, 3*b+1:4*b) = rand(b,b);
%%
% sanity check: 
rank(A)
rank(A(1:2*b, 2*b+1:end))
rank(A(2*b+1:end, 1:2*b))
rank(A(1:b, b+1:2*b))
rank(A(b+1:2*b, 1:b))
rank(A(2*b+1:3*b, 3*b+1:4*b))
rank(A([1:8, 17:32], 9:16))
rank(A(1:b,1:b))
%%
% now construct the HSS approximation: 
H = hss(A,blocksize = 8,k = r);
%%
% testing blockbuilder
A12 = blockbuilder(H.A12,H);
trueA12 = A(1:32,33:64);
A12-trueA12
%%
% testing subsef
% only in A11
norm(H(5:14,3:7) - A(5:14,3:7))/norm(A(5:14,3:7))
% only in A22
norm(H(45:54,43:57) - A(45:54,43:57))/norm(A(45:54,43:57))
% only in A21
norm(H(45:54,2:10) - A(45:54,2:10))/norm(A(45:54,2:10))
% only in A12
norm(H(17:26,43:57) - A(17:26,43:57))/norm(A(17:26,43:57))

% A11 and A12
norm(H(17:26,23:57) - A(17:26,23:57))/norm(A(17:26,23:57))
% A11 and A21
norm(H(22:54,2:10) - A(22:54,2:10))/norm(A(22:54,2:10))
% A12 and A22
norm(H(13:54,43:57) - A(13:54,43:57))/norm(A(13:54,43:57))
% A21 and A22
norm(H(45:54,2:60) - A(45:54,2:60))/norm(A(45:54,2:60))

% all 4
norm(H(30:38,11:60) - A(30:38,11:60))/norm(A(30:38,11:60))


%% test off diag leaves
idxr = 1:b; 
idxc = b+1:2*b; 
L = leafbuild(H.A11.A12);
norm(L-A(idxr, idxc))
%%
idxc = 1:b; 
idxr = b+1:2*b; 
L = leafbuild(H.A11.A21);
norm(L-A(idxr, idxc))
%%
idxr = 2*b + 1: 3*b; 
idxc = 3*b+1:4*b; 
L = leafbuild(H.A22.A12);
norm(L-A(idxr, idxc))
%%
idxc = 2*b + 1: 3*b; 
idxr = 3*b+1:4*b; 
L = leafbuild(H.A22.A21);
norm(L-A(idxr, idxc))
%%
% check diags: 
LL = {H.A11.A11, H.A11.A22, H.A22.A11, H.A22.A22};
for k = 1:4
    L = LL{k}; 
    L = leafbuild(L);
    idx = (k-1)*b+1:k*b; 
    err(k) = norm(L-A(idx,idx))
end
%%
% now check the next level up
norm(blkdiag(H.A11.A12.Z,H.A11.A21.Z)*H.A12.Z*H.A12.lrcomponent*H.A12.Y*blkdiag(H.A22.A21.Y,H.A22.A12.Y)-A(1:16,17:32))
norm(blkdiag(H.A22.A12.Z,H.A22.A21.Z)*H.A21.Z*H.A21.lrcomponent*H.A21.Y*blkdiag(H.A11.A21.Y,H.A11.A12.Y)-A(17:32,1:16))

%%
% check matvec: 
v = rand(m,1);
mult_err = norm(H*v-A*v)









