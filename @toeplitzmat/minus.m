function h = minus(T, G)
%-  Minus for TOEPLITZMAT objects
%  T - G subtracts G from T. 
%  T and G can be scalars or matrices or TOEPLITZMATs

ftest = @(x) (isnumeric(x) || islogical(x));
h = toep_apply(@(x,y) x - y, ftest, "-", T, G);
end

