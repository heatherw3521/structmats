function h = lt(T,G)
%<  Pointwise less than for DIAGONALSTRUCTURE objects
%   T and g can be DIAGONALSTRUCTUREs or scalars
%   behaves the same as < for matrices in MATLAB

ftest = @(x) (isnumeric(x) || islogical(x));
h = diag_apply(@(x,y) x < y, ftest, "<", T, G);
end