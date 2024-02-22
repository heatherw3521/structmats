function h = ldivide(G, T)
%.\   Pointwise TOEPLITZMAT left divide
%   G .\ T pointwise divides T by g. 
%   T and g can be scalars or Matrices or TOEPLITZMATs

h = rdivide(T,G);
end

