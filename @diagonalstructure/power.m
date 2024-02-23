function h = power(T, G)
%.^  Pointwise power for DIAGONALSTRUCTURE objects
%  T .^ G raies T to the G, pointwise. 
%  T and G can be scalars or Matrices or DIAGONALSTRUCTUREs

ftest = @(x) (isnumeric(x) || islogical(x));
h = diag_apply(@(x,y) x .^ y, ftest, ".^", T, G);
end

