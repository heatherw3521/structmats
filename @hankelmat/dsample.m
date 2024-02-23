function h = dsample(H,i,j)
%SAMPLE Summary of this function goes here
%   Detailed explanation goes here

n = size(H,1);
h = (i+j<=n) .* H.hc(i+j-1).' + (i+j>n).*H.hr(i+j-n);

end

