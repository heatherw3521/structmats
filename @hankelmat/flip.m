function h = flip(H)
%FLIP up-down flip for HANKELMAT
%   The flip of a Hankel matrix is a Toeplitz matrix!

h = toeplitzmat(flip(H.hc), H.hr);
end

