function h = ldivide(G, T)
%.\   Pointwise DIAGONALSTRUCTURE left divide
%   G .\ T pointwise divides T by g. 
%   T and g can be scalars or Matrices or DIAGONALSTRUCTUREs

h = rdivide(T,G);
end

