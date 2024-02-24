function h = dsample(H,i,j)
%SAMPLE Summary of this function goes here
%   Detailed explanation goes here

m = size(H,1); n = size(H,2); k = max(n,m);
h = (i+j-1<=m).*H.hc(min(i+j-1, m)).' + (i+j-1>m).*H.hr(min(abs(i+j-m-1)+1,n));

end

