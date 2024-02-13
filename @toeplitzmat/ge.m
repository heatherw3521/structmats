function h = ge(T,g)
%>=  Pointwise greater than or equal to for TOEPLITZMAT objects
%   T and g can be TOEPLITZMATs or scalars
%   behaves the same as >= for matrices in MATLAB

ftest = @(x) isnumeric(x);
h = toepcompare(T, g, @(x,y) x >= y, ">=", ftest);
end