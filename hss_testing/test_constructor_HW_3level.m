% test constructor on an example for a 3 level tree
clear all
close all
%%
warning off
m = 2^6; n = 2^6;
blocksize = 8;
r = 5; 

X = linspace(0,1,n);
Y = linspace(4,5,m).';
C = randn(m, r)*randn(r,n);
C = C + eye(m,n); 
%% 
% leaf
H = hss(C,blocksize = blocksize,k = r);
L = leafbuild(H.A11.A11.A12);
idxr = 1:8; 
idxc = 9:16; 
norm(L-C(idxr, idxc))
%%
%leaf
H = hss(C,blocksize = blocksize,k = r);
L = leafbuild(H.A22.A22.A12);
idxr = 49:56; 
idxc = 57:64; 
norm(L-C(idxr, idxc))
%%
% leaf
H = hss(C,blocksize = blocksize,k = r);
L = leafbuild(H.A22.A11.A21);
idxr = 41:48;
idxc = 33:40; 
norm(L-C(idxr, idxc))
%%
% level 2 block
norm(blkdiag(H.A11.A11.A12.Z,H.A11.A11.A21.Z)*H.A11.A12.Z*H.A11.A12.lrcomponent...
    *H.A11.A12.Y*blkdiag(H.A11.A22.A21.Y,H.A11.A22.A12.Y)-C(1:16,17:32))

%%
% level 2 block
norm(blkdiag(H.A22.A22.A12.Z,H.A22.A22.A21.Z)*H.A22.A21.Z*H.A22.A21.lrcomponent...
    *H.A22.A21.Y*blkdiag(H.A22.A11.A21.Y,H.A22.A11.A12.Y)-C(49:64,33:48))
%%
% level 3 block: 


