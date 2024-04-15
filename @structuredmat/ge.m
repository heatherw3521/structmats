function h = ge(A,B)
%>=  Pointwise greater than or equal to for STRUCTUREDMAT objects
%   T and g can be STRUCTUREDMATs or scalars
%   behaves the same as >= for matrices in MATLAB

warning(['Using Naive implementation of >= for inputs of type %s and %s.', ...
           'This should be overloaded to take advantage of the structures!'], ...
           class(A), class(B));
h = ge(full(A), full(B));
end