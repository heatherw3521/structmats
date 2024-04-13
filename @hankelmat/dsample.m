function h = dsample(H,i,j)
%SAMPLE Samples from a hankel matrix, given subindexing sets, i,j. This is
%   done efficiently without forming the full matrix.

%% TODO: test this, there might be an issue like in dsample for toeplitz...

m = size(H,1); n = size(H,2); k = max(n,m);
h = (i+j-1<=m).*H.hc(min(i+j-1, m)).' + (i+j-1>m).*H.hr(min(abs(i+j-m-1)+1,n));

end

