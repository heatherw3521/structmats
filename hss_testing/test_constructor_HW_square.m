%test constructor and matvec for square powers of 2: 
clear all
close all
m = 2^10; n = 2^10; 
blocksize = 8; 
r = 4; 
C = randn(m,r)*randn(r,n);
C = C + eye(m);

%%
H = hss(C,blocksize = blocksize,k = r);
%%
v = randn(n,1);
norm(H*v-C*v)

%%
% test constructor and matvec for non powers of 2. Still square
clear all
m = 2^10+5; n = 2^10+5; 
blocksize = 8; 
r = 4; 
C = randn(m,r)*randn(r,n);
C = C + eye(m);

%%
H = hss(C,blocksize = blocksize,k = r);
%%
v = randn(n,1);
norm(H*v-C*v)