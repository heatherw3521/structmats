function Y = mtimes(H,X)
%MTIMES Matrix Times for HANKELMAT objects
%   H,X can be Matrices, Scalars, or Hankelmat's

% We compute the product efficiently by noting that the up-down flip
%   of a Hankel Matrix is Toeplitz
    
% TODO: Maybe just hard-code this to avoid toeplitzmat instantiation
% overhead
if(isa(H,"hankelmat"))
    Y = flip(toeplitzmat(flip(H.hc),H.hr) * X);
else
    Y = (X.' * H.').';
end
end

