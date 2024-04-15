function h = times(A, B)
%.*   Pointwise multiplication for STRUCTUREDMAT objects
%   A .* B pointwise multiplies A and B. 
%   A and B can be scalars or Matrices or STRUCTUREDMATs

warning(['Using Naive implementation of .* for inputs of type %s and %s.', ...
           'This should be overloaded to take advantage of the structures!'], ...
           class(A), class(B));
h = times(full(A),full(B));
end