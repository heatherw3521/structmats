function h = times(T, G)
%.*   Pointwise multiplication for TOEPLITZMAT objects
%   T .* G pointwise multiplies T and G. 
%   T and G can be scalars or Matrices or TOEPLITZMATs

ftest = @(x) (isnumeric(x) || islogical(x));
h = toep_apply(@(x,y) x .* y, ftest, ".*", T, G);
end