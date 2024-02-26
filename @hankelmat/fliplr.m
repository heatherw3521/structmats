function h = fliplr(H)
%FLIP left-right flip for HANKELMAT
%   The flip of a Hankel matrix is a Toeplitz matrix!

Hflat = [H.hc.' H.hr(2:end)];
tc = Hflat(size(H,2):end);
tr = flip(Hflat(1:size(H,2)));
h = toeplitzmat(tc, tr);
end

