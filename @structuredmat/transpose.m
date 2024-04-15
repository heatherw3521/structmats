function v = ctranspose(A)
%CTRANSPOSE Conjugate Transpose for STRUCTUREDMAT objects

warning(['Using Naive implementation of TRANSPOSE for input of type %s.', ...
           'This should be overloaded to take advantage of the structure!'], ...
           class(A));
v = transpose(full(A));
end

