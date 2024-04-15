function h = mtimes(A, B)
%MTIMES Matrix Times for STRUCTUREDMAT objects

warning(['Using Naive implementation of * for inputs of type %s and %s.', ...
           'This should be overloaded to take advantage of the structures!'], ...
           class(A), class(B));
h = mtimes(full(A),full(B));
end

