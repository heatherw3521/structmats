%Examples_Solvers

%This package will become a collection of fast solvers for highly
%structured matrices. These currently include: 
%
% 1) Rectangular vandermonde matrices with nodes on the unit circle (Non-uniform
% discrete Fourier transform of type-II). 
%
% 2) Square Toeplitz matrices.
%%

%% Example 1: Nonuniform discrete Fourier transform
%
% We solve the system Vx = b, where V is an m by n Vandermonde matrix with 
% entries V(j,k) = exp(-2*pi*1i*p(j)*(k-1)). The m by 1 vector p contains 
% sample locations on [0, 1]. These do not have to be ordered. 
%
% We first set up the problem: 
n = 2050; % number of modes
m = 2*n; 
mds = 0:n-1;
p = rand(m,1); %sample locations in [0, 1] for non-unif DFT 
V = exp(-2*pi*1i*p.*mds); b = V*rand(n,1); %get a RHS in range of V. 
%%
% To solve with structmat, we only need the node locations, RHS, and number
% of modes. Here we call the solver: 

x = structsolv_nudft2(p, n, b); 
%%
% check residual: 
norm(V*x-b)./norm(b)
%
% More nonuniform DFT examples, including with multiple RHSs and varying 
% error tolerances, are given in Example_NUDFT.m
%
%% Example 2: Toeplitz system
%
% A Toeplitz matrix is a matrix with constant diagonal entries. MATLAB has
% a built in Toeplitz command we can use to build a sample Toeplitz matrix.
% It requires the first row (tr), and first column (tc) of the matrix. 

n = 2040; %size of Toeplitz matrix       
tr = rand(n,1); 
tc = rand(n,1);  
x = rand(n,1); %true solution
b = T*x; %RHS

%%
% We can solve the system Tx = b in near-linear time with our fast Toeplitz
% solver. It only requires tr, tc, and the right hand side. 

xns = structsolv_toeplitz(tc, tr,b,'tol',tol );
%%
% check error
norm(x-xns)/norm(x)
%
% For more examples with Toeplitz matrices, see Example_Toeplitz.m




