function T = nearby_toeplitz(A)
%NEARBY_TOEPLITZ Finds the closest Toeplitz matrix to A, in the sense of
%   the Frobenius norm. It can be shown that the closest matrix is obtained
%   by simply averaging along the diagonals!
%
%   TODO: Add proof to devlog. It's just basic calculus.

[m,n] = size(A);

% Compute the first column
tc = zeros(1,m);
for i = 1:m
    tc(i) = mean(diag(A,1-i));
end

% Compute the first row
tr = zeros(1,n);
tr(1) = tc(1);
for j = 2:n
    tr(j) = mean(diag(A,j-1));
end

T = toeplitzmat(tc,tr);

end

