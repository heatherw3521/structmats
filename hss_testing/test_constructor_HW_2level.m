%square test: this is  a simple test to see if the constructor does 
% the right thing for a two-level tree. 

clear all
close all
%%
m = 47; 
n = 47; 
bl = 15; % blocksize
r = 4; % rank
% x = linspace(0,1,m);
% y = x+0.01;
% A = gallery('cauchy',x,-y);
% A = A(1:m,1:n);
global A
A = rand(m,r)*rand(r,n) + eye(m,n);
%A = rand(m,n);
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
rank(A(1:2*bl, 2*bl+1:end))
rank(A(2*bl+1:end, 1:2*bl))
rank(A(1:bl, bl+1:2*bl))
rank(A(bl+1:2*bl, 1:bl))
rank(A(2*bl+1:3*bl, 3*bl+1:4*bl))
rank(A([1:8, 17:32], 9:16))
rank(A(1:bl,1:bl))
%%
bl = 10;
r = 4;
for m = 82:200
    for n = 50:200
        A = rand(m,r)*rand(r,n);
        H = hss(A,blocksize = bl,k = r);
    end
end
%%
% now construct the HSS approximation: 
H = hss(A,blocksize = bl);
b = ones(m,1);
global x
x = A\b;
%% 
% testing constructor on weird sizes
spy(H)
%%
% testing blockbuilder
A12 = blockbuilder(H.A12,H);
trueA12 = A(1:ceil(m/2),ceil(n/2)+1:end);
A12-trueA12
norm(A12-trueA12)
%%
% blockbuilder lower level
block = H.A11.A12
estfullblock = blockbuilder(block,H.A11)
truefullblock = (A(block.Ir(1):block.Ir(2),block.Ic(1):block.Ic(2)));
size(estfullblock)
size(truefullblock)
estfullblock-truefullblock
norm(estfullblock-truefullblock)/norm(truefullblock)
%%
% one more level down
block = H.A22.A11.A21
estfullblock = blockbuilder(block,H.A22.A11)
truefullblock = (A(block.Ir(1):block.Ir(2),block.Ic(1):block.Ic(2)));
size(estfullblock)
size(truefullblock)
estfullblock-truefullblock
norm(estfullblock-truefullblock)/norm(truefullblock)
%%
% working with just the block above (leaf level) 
block.Z*block.lrcomponent*block.Y - truefullblock
[Z,rows] = inter_decomp_tfile(truefullblock,orientation = 'rows')
[Y,cols] = inter_decomp_tfile(truefullblock,orientation = 'columns')
lrcomponent = truefullblock(rows,cols)
Z*lrcomponent*Y - truefullblock
%%
% above works lets try decomping the whole row and col
[Z,rows] = inter_decomp_tfile(A(block.Ir(1):block.Ir(2),:),orientation = 'rows')
[Y,cols] = inter_decomp_tfile(A(:,block.Ic(1):block.Ic(2)),orientation = 'columns')
lrcomponent = truefullblock(rows,cols)
Z*lrcomponent*Y - truefullblock



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
idxr = 1:bl; 
idxc = bl+1:2*bl; 
L = leafbuild(H.A11.A12);
norm(L-A(idxr, idxc))
%%
idxc = 1:bl; 
idxr = bl+1:2*bl; 
L = leafbuild(H.A11.A21);
norm(L-A(idxr, idxc))
%%
idxr = 2*bl + 1: 3*bl; 
idxc = 3*bl+1:4*bl; 
L = leafbuild(H.A22.A12);
norm(L-A(idxr, idxc))
%%
idxc = 2*bl + 1: 3*bl; 
idxr = 3*bl+1:4*bl; 
L = leafbuild(H.A22.A21);
norm(L-A(idxr, idxc))
%%
% check diags: 
LL = {H.A11.A11, H.A11.A22, H.A22.A11, H.A22.A22};
for k = 1:4
    L = LL{k}; 
    L = leafbuild(L);
    idx = (k-1)*bl+1:k*bl; 
    err(k) = norm(L-A(idx,idx))
end
%%
% now check the next level up
norm(blkdiag(H.A11.A12.Z,H.A11.A21.Z)*H.A12.Z*H.A12.lrcomponent*H.A12.Y*blkdiag(H.A22.A21.Y,H.A22.A12.Y)-A(1:16,17:32))
norm(blkdiag(H.A22.A12.Z,H.A22.A21.Z)*H.A21.Z*H.A21.lrcomponent*H.A21.Y*blkdiag(H.A11.A21.Y,H.A11.A12.Y)-A(17:32,1:16))

%%
% check matvec: 
v = rand(n,1);
mult_err = norm(H*v-A*v)









