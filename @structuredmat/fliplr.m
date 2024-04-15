function v = fliplr(A)
%CTRANSPOSE Left-right flip for STRUCTUREDMAT objects

warning(['Using Naive implementation of FLIPLR for input of type %s.', ...
           'This should be overloaded to take advantage of the structure!'], ...
           class(A));
v = flip(full(A));
end