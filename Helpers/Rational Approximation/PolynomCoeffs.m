function c = PolynomCoeffs(Z)
%POLYNOMCOEFFS Gets the coefficients of a monomial with roots given by Z.
%   Does this by iteratively building up the coefficients. The complexity
%   is O(n^2), and I'm not sure if it can be improved.
arguments
    Z (:,1) {mustBeNumeric}
end

n = length(Z);
c = [1 ; zeros([n 1])];
for d = 1:n
    c = -Z(d)*c + [0 ; c(1:end-1)];
end

end

