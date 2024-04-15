function h = and(A,B)
%&  Pointwise logical and for STRUCTUREDMAT objects
%   A and B can be STRUCTUREDMATs or scalars
%   behaves the same as & for matrices in MATLAB

warning(['Using Naive implementation of & for inputs of type %s and %s.', ...
           'This should be overloaded to take advantage of the structures!'], ...
           class(A), class(B));
h = and(full(A), full(B));
end