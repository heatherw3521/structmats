function C = nearby_circulant(A)
%NEARBY_CIRUCLANT Finds the closest Circulant matrix to A, in the sense of
%   the Frobenius norm. It can be shown that the closest matrix is obtained
%   by simply averaging along the diagonals corresponding a constant
%   value in the Circulant matrix. For instance, the second element of the
%   first column will be the total average of the first subdiagonal and
%   last superdiagonal.
%
%   TODO: Add proof to devlog. It's just basic calculus.

[m,n] = size(A);
if(m~=n)
    error(['Cannot project a non-square (%s-by-%s) matrix to ', ...
        'a circulant matrix!'], ""+m,""+n);
end

% Compute the first column
tc = zeros(1,m);
tc(1) = mean(diag(A,0));
for i = 2:m
    tc(i) = (sum(diag(A,1-i))+sum(diag(A,m-i+1)))/m;
end
C = circulantmat(tc);

end

