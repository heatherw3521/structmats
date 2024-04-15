function h = uminus(A)
%-A  Unary Minus for STRUCTUREDMAT objects
%   -A negates A pointwise, where A is a STRUCTUREDMAT object

warning(['Using Naive implementation of uminus (-A) for input of type %s.', ...
           'This should be overloaded to take advantage of the structures!'], ...
           class(A));
h = -1 .* A;
end

