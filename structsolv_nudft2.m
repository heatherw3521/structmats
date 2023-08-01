%structmat solver: 
function varargout = structsolv_nudft2(p,N,b, varargin)
% solve a least squares problem Vx = b, where V is 1D type-II nonuniform DFT matrix,
% i.e., V(j,k) = exp(-2i*pi*p(j)*(k-1))
% p = locations on [0, 1] where signal is sampled.
% N = number of modes. (modes are 0:N-1)
% b = right-hand side (can also be a matrix of multiple right hand sides)
%
% x = structsolv_nudft2(nodes, n, b, 'tol', tol) solves to the tolerance tol.
% default output is x, which is an approximate least-squares solution.
%
% [L, P, x] = structsolv_nudft2(nodes, ...) returns L, a urv factorization object for a permutation
% of V. P is the permutation information, and x is the solution. L and P
% can be used in apply_inudft_solve(...) to solve Vx = y for additional right-hand sides. 
% 
% Our method is based on [1]. 

% References: 
% [1] Barnett, Epperly, Wilber, "A superfast direct inversion method for the nonuniform
% discrete Fourier transform" (not yet public)

%%
% convert to Vandermonde nodes, and then call NUDFT code
ucnodes = exp(-2i*pi*p(:));
varargout = INUDFT(ucnodes,N, b, varargin);
varargout = varargout{:}; 
end


