function h = mpower(T,G)
%MPOWER Matrix Power for TOEPLITZMAT objects.
%   T and G can be TOEPLITZMATs, scalars, or matrices
%
%   I'm not sure if there's an efficient way to do this,
%   so the code below just casts to regular matrices.

% Attempt matrix exponential
h = full(T) ^ full(G);
end

