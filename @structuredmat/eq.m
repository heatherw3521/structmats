function h = eq(A,B)
%==  Pointwise equal to for STRUCTUREDMAT objects
%   A and B can be STRUCTUREDMATs or scalars
%   behaves the same as == for matrices in MATLAB

warning(['Using Naive implementation of == for inputs of type %s and %s.', ...
           'This should be overloaded to take advantage of the structures!'], ...
           class(A), class(B));
eq(full(A),full(B));
end