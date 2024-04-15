function h = ldivide(A, B)
%.\   Pointwise STRUCTUREDMAT left divide
%   A .\ B pointwise divides A by B. 
%   A and B can be scalars or Matrices or STRUCTUREDMATs

% Simply reuse rdivide, with the proper call
%   i.e.    A .\ B = B ./ A
warning(['Using Naive implementation of .\ for inputs of type %s and %s.', ...
           'This should be overloaded to take advantage of the structures!'], ...
           class(A), class(B));
h = rdivide(B,A);
end

