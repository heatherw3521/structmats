function h = ldivide(g, T)
%./   Pointwise TOEPLITZMAT left divide
%   g ./ T pointwise divides T by g. 
%   T and g can be scalars or TOEPLITZMATs

h = rdivide(T,g);
end

