function h = power(A, B)
%.^  Pointwise power for STRUCTUREDMAT objects
%  A .^ B raies A to the B, pointwise. 
%  A and B can be scalars or Matrices or STRUCTUREDMATs

warning(['Using Naive implementation of .^ for inputs of type %s and %s.', ...
           'This should be overloaded to take advantage of the structures!'], ...
           class(A), class(B));
h = power(full(A),full(B));
end

