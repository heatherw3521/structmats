function h = mrdivide(T, G)
%MRDIVIDE Matrix Right Divide for TOEPLITZMAT objects.
%   T and G can be TOEPLITZMATs, scalars, or matrices
%
%   The code below is just a placeholder!

% Basic implementation that does not take advantage of structure.
% The way that this will be efficiently overloaded depends on the type
% of diagonally structured matrix
h = full(T) / full(G);
end

