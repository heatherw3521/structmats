function h = not(A)
%~  Pointwise logical not for STRUCTUREDMAT objects
%   A is a STRUCTUREDMAT
%   behaves the same as ~ for matrices in MATLAB

warning(['Using Naive implementation of ~ for input of type %s.', ...
           'This should be overloaded to take advantage of the structure!'], ...
           class(A));
h = not(full(A));
end