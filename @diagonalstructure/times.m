function h = times(T, G)
%.*   Pointwise multiplication for DIAGONALSTRUCTURE objects
%   T .* G pointwise multiplies T and G. 
%   T and G can be scalars or Matrices or DIAGONALSTRUCTUREs

ftest = @(x) (isnumeric(x) || islogical(x));
h = diag_apply(@(x,y) x .* y, ftest, ".*", T, G);
end