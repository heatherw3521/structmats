function h = dsample(T,i,j)
%SAMPLE Summary of this function goes here
%   Detailed explanation goes here

h = (j>=i) .* T.tr(min(abs(j-i)+1,T.size(2))) + (j<i).*T.tc(min(abs(i-j)+1,T.size(1))).';

end

