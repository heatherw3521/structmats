function X = Toep_solve(L,B)
%solves TX = B where L stores urv factorization for the transformed 
% Toeplitz matrix T.  B is multiple RHS. 
%
% Example: First call [L, ~] = structsolv_toeplitz(tr,tc, n, b1); to 
% produce L, where tr, tc are the first row, col (resp) of T.  
% Now solve system TX = B with X = Toep_solve(L,p, B); 


% TO DO: Optimize urv_solve wrt BLAS3 ops
[n,s] = size(B); 
B = sqrt(n)*ifft(B); 
X = zeros(L.n, s);
for j= 1:s
    X(:,j) = urv_solve(L, B(:,j));
end
N = (1:n).'; 
w = exp(pi*1i/n);
D0 = spdiags(w.^(-(N-1)), 0, n,n);
X = D0*fft(X)/sqrt(n);

end
