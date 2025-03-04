function AA = buildcauchy_old(J, K, N, nodes,G, L)
%build the cauchy matrix from the given indices:
J = J(:); 
K = K(:); 
w = exp(1i*pi/N); 
%a = w.^(2*(J-1));
a = nodes(J); 
%b = w.^(2*K-1); 
b = w.^(2*K); b = b(:);
A = bsxfun(@minus, a, b.'); 
A = 1./A; 

%apply hadamard product
[~,r] =size(G); 
k = length(J);
n = length(K); 
AA = 0; 
for j = 1:r
    B = spdiags(G(J, j), 0,k,k)*A*spdiags(conj(L(K,j)), 0, n,n);
    AA = AA +B; 
end