function h = mtimes(T, G)
%MTIMES Matrix Times for DIAGONALSTRUCTURE objects
%   This is where the magic happens! Matrix-matrix products can be
%   computed much cheaper when one of the two is a toeplitz matrix.

% Basic implementation that does not take advantage of structure.
% The way that this will be efficiently overloaded depends on the type
% of diagonally structured matrix
h = full(T) * full(G);
end

