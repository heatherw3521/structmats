function h = and(T,G)
%&  Pointwise logical and for DIAGONALSTRUCTURE objects
%   T and g can be DIAGONALSTRUCTUREs or scalars
%   behaves the same as & for matrices in MATLAB

% we can do logic on numerical data too
ftest = @(x) (isnumeric(x) || islogical(x));
h = diag_apply(@(x,y) x & y, ftest, "&", T, G);
end