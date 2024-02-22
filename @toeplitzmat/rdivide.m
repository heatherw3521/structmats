function h = rdivide(T, G)
%./   Pointwise TOEPLITZMAT right divide
%   T ./ G pointwise divides T by G. 
%   T and G can be scalars or matrices or TOEPLITZMATs


ftest = @(x) (isnumeric(x) || islogical(x));
h = toep_apply(@(x,y) x ./ y, ftest, "./", T, G);
end

