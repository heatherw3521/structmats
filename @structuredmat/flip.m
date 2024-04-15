function v = flip(A)
%CTRANSPOSE Up-down flip for STRUCTUREDMAT objects

warning(['Using Naive implementation of FLIP for input of type %s.', ...
           'This should be overloaded to take advantage of the structure!'], ...
           class(A));
v = flip(full(A));
end