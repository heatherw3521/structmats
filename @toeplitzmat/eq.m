function h = eq(T,g)
%==  Pointwise equal to for TOEPLITZMAT objects
%   T and g can be TOEPLITZMATs or scalars
%   behaves the same as == for matrices in MATLAB

h = toepcompare(T, g, @(x,y) x == y, "==");
end