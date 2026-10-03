function [b, xs] = mp_manufactured(H, iscomplex)
%MP_MANUFACTURED  Right-hand side with a KNOWN minimum-norm solution:
%   z ~ N(0,I),  xs = H' z  (in Range(H^*)),  b = H xs.
%   Then xs is exactly the minimum-norm solution of H x = b, so
%   ||x - xs||/||xs|| measures the solver alone -- with O(N) HSS products,
%   no dense matrix, at any size.
if nargin < 2, iscomplex = false; end
m = H.size(1);
z = randn(m,1);
if iscomplex, z = z + 1i*randn(m,1); end
xs = H' * z;
b = H * xs;
end
