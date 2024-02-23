function h = mldivide(T, G)
%MLDIVIDE Matrix Left Divide for TOEPLITZMAT objects.
%   T and G can be TOEPLITZMATs, scalars, or matrices
%
%   The code below is just a placeholder!

if(isa(T,'toeplitzmat'))
    T = toeplitz(T);
end
if(isa(G,'toeplitzmat'))
    G = toeplitz(G);
end
% Attempt matrix division
h = T \ G;
end

