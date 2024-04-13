function h = mldivide(T, G)
%MLDIVIDE Matrix Left Divide for DIAGONALSTRUCTURE objects.
%   T and G can be DIAGONALSTRUCTUREs, scalars, or matrices
%
%   The code below is just a placeholder!

% Basic implementation that does not take advantage of structure.
% The way that this will be efficiently overloaded depends on the type
% of diagonally structured matrix
h = full(T) \ full(G);
end

