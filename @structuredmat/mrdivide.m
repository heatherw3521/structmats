function h = mrdivide(A, B)
%MRDIVIDE Matrix Right Divide for STRUCTUREMAT objects.
%   A and B can be STRUCTUREDMATs, scalars, or matrices

warning(['Using Naive implementation of / for inputs of type %s and %s.', ...
           'This should be overloaded to take advantage of the structures!'], ...
           class(A), class(B));
h = mrdivide(full(A),full(B));
end

