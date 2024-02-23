function h = mtimes(T, G)
%MTIMES Matrix Times for DIAGONALSTRUCTURE objects
%   This is where the magic happens! Matrix-matrix products can be
%   computed much cheaper when one of the two is a toeplitz matrix.

% Attempt the multiplication in an inefficient but consistent way.
% This should be overloaded!
h = full(T) * full(G);
end

