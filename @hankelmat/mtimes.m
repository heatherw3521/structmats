function Y = mtimes(H,X)
%MTIMES Matrix Times for HANKELMAT objects
%   H,X can be Matrices, Scalars, or Hankelmat's

% We compute the product efficiently by noting that the up-down flip
%   of a Hankel Matrix is Toeplitz
    
% TODO: Maybe just hard-code this to avoid toeplitzmat instantiation
% overhead
if(isa(H,"hankelmat"))

    if(isequal(H.TH, "uncomputed"))
        H.TH = toeplitzmat(flip(H.hc),H.hr);
    end
    
    Y = flip(H.TH * X);
else
    Y = (X.' * H.').';
end
end

