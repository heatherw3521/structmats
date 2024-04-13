function h = rdivide(T, G)
%./   Pointwise DIAGONALSTRUCTURE right divide
%   T ./ G pointwise divides T by G. 
%   T and G can be scalars or matrices or DIAGONALSTRUCTUREs

ftest = @(x) (isnumeric(x) || islogical(x));
h = diag_apply(@(x,y) x ./ y, ftest, "./", T, G);
end

