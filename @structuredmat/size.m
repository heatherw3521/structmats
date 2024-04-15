function v = size(A)
%SIZE Size for STRUCTUREDMAT objects

warning(['Using Naive implementation of SIZE for input of type %s.', ...
           'This should be overloaded to take advantage of the structure!'], ...
           class(A));
v = size(full(A));

