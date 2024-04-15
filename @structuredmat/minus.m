function h = minus(A, B)
%-  Minus for STRUCTUREDMAT objects
%  A - B subtracts B from A. 
%  A and B can be scalars or matrices or STRUCTUREDMATs

warning(['Using Naive implementation of > for inputs of type %s and %s.', ...
           'This should be overloaded to take advantage of the structures!'], ...
           class(A), class(B));
h = minus(full(A),full(B));
end

