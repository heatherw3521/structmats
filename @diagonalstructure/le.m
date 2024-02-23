function h = le(T,G)
%<=  Pointwise less than or equal to for DIAGONALSTRUCTURE objects
%   T and g can be DIAGONALSTRUCTUREs or scalars
%   behaves the same as <= for matrices in MATLAB

ftest = @(x) (isnumeric(x) || islogical(x));
h = diag_apply(@(x,y) x <= y, ftest, "<=", T, G);
end