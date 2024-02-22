function h = le(T,G)
%<=  Pointwise less than or equal to for TOEPLITZMAT objects
%   T and g can be TOEPLITZMATs or scalars
%   behaves the same as <= for matrices in MATLAB

ftest = @(x) (isnumeric(x) || islogical(x));
h = toep_apply(@(x,y) x <= y, ftest, "<=", T, G);
end