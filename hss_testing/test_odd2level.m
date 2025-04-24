%square test: this is  a simple test to see if the constructor does 
% the right thing for a two-level tree. 

clear all
close all
%%
m = 37; 
n = m; 
b = 10; % blocksize
% if b = 9 we see issues with the matvec. if b = 10 this is corrected.
% when b = 9 there is a single level H.A11.A11 where the block is 10 by 10,
% so it gets cut into a 5 by 5 collection (level 3). However, the other
% level 2 blocks are not cut smaller so this causes some issues in the
% matvec. It is avoided if one uses 


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
H = hss(A,blocksize = b,k = r);
%% test off diag leaves
idxr = 1:10; 
idxc = 11:19; 
L = leafbuild(H.A11.A12);
norm(L-A(idxr, idxc))
%%
idxc = 1:10; 
idxr = 11:19; 
L = leafbuild(H.A11.A21);
norm(L-A(idxr, idxc))
%%
idxr = 20:28; 
idxc = 29:37; 
L = leafbuild(H.A22.A12);
norm(L-A(idxr, idxc))
%%
idxc = 20:28;
idxr = 29:37; 
L = leafbuild(H.A22.A21);
norm(L-A(idxr, idxc))
%%
%%
% now check the next level up
norm(blkdiag(H.A11.A12.Z,H.A11.A21.Z)*H.A12.Z*H.A12.lrcomponent*H.A12.Y*blkdiag(H.A22.A21.Y,H.A22.A12.Y)-A(1:19,20:37))
norm(blkdiag(H.A22.A12.Z,H.A22.A21.Z)*H.A21.Z*H.A21.lrcomponent*H.A21.Y*blkdiag(H.A11.A21.Y,H.A11.A12.Y)-A(20:37,1:19))

%%
% check matvec: 
v = rand(m,1);
mult_err = norm(H*v-A*v)









