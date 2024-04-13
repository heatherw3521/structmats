function h = fliplr(T)
%FLIP left-right flip for TOEPLITZMAT
%   The flip of a Toeplitz matrix is a Hankel matrix!

Tflat = [flip(T.tr(2:end)) T.tc.'];
hc = Tflat(1:size(T,1)); 
hr = Tflat(size(T,1):end);
h = hankelmat(hc, hr);
end

