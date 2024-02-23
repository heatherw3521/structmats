function h = minus(T, G)
%-  Minus for DIAGONALSTRUCTURE objects
%  T - G subtracts G from T. 
%  T and G can be scalars or matrices or DIAGONALSTRUCTUREs

ftest = @(x) (isnumeric(x) || islogical(x));
h = diag_apply(@(x,y) x - y, ftest, "-", T, G);
end

