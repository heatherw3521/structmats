function h = gt(T,G)
%>  Pointwise greater than for DIAGONALSTRUCTURE objects
%   T and g can be DIAGONALSTRUCTUREs or scalars
%   behaves the same as > for matrices in MATLAB

ftest = @(x) (isnumeric(x) || islogical(x));
h = diag_apply(@(x,y) x > y, ftest, ">", T, G);
end