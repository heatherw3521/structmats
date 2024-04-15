function v = ctranspose(A)
%CTRANSPOSE Conjugate Transpose for STRUCTUREDMAT objects

warning(['Using Naive implementation of CTRANSPOSE for input of type %s.', ...
           'This should be overloaded to take advantage of the structure!'], ...
           class(A));
v = ctranspose(full(A));
end

