function h = ne(T,g)
%~=  Pointwise not equal to for TOEPLITZMAT objects
%   T and g can be TOEPLITZMATs or scalars
%   behaves the same as ~= for matrices in MATLAB

h = toepcompare(T, g, @(x,y) x ~= y, "~=");
end