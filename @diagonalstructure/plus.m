function h = plus(T, G)
%+   Plus for DIAGONALSTRUCTURE objects
%  T + G adds T and G. 
%  T and G can be scalars or Matrices or DIAGONALSTRUCTUREs

ftest = @(x) (isnumeric(x) || islogical(x));
h = diag_apply(@(x,y) x + y, ftest, "+", T, G);
end

