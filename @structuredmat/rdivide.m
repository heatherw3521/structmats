function h = rdivide(A, B)
%./   Pointwise STRUCTUREDMAT right divide
%   A ./ B pointwise divides A by B. 
%   A and B can be scalars or matrices or STRUCTUREDMATs

warning(['Using Naive implementation of ./ for inputs of type %s and %s.', ...
           'This should be overloaded to take advantage of the structures!'], ...
           class(A), class(B));
h = rdivide(full(A),full(B));
end

