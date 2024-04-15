function h = mpower(A,B)
%MPOWER Matrix Power for STRUCTUREDMAT objects.
%   A and B can be STRUCTUREDMATs, scalars, or matrices

warning(['Using Naive implementation of ^ for inputs of type %s and %s.', ...
           'This should be overloaded to take advantage of the structures!'], ...
           class(A), class(B));
h = mpower(full(A),full(B));
end

