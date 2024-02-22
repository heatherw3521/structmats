function h = power(T, G)
%.^  Pointwise power for TOEPLITZMAT objects
%  T .^ G raies T to the G, pointwise. 
%  T and G can be scalars or Matrices or TOEPLITZMATs

ftest = @(x) (isnumeric(x) || islogical(x));
h = toep_apply(@(x,y) x .^ y, ftest, ".&", T, G);
end

