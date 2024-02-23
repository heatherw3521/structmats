function h = and(T,G)
%|  Pointwise logical and for TOEPLITZMAT objects
%   T and g can be TOEPLITZMATs or scalars
%   behaves the same as & for matrices in MATLAB

% we can do logic on numerical data too
ftest = @(x) (isnumeric(x) || islogical(x));
h = toep_apply(@(x,y) x & y, ftest, "&", T, G);
end