function h = dsample(T,I,J)
%SAMPLE Summary of this function goes here
%   Detailed explanation goes here

C = T.tc(min(abs(I-J)+1,size(T,1)));
R = T.tr(min(abs(J-I)+1,size(T,2)));

% MATLAB doesn't retain dimensions, in some cases... great.
if(length(J) == 1) R = R.'; end
if(length(I) == 1) C = C.'; end
% This whole function would be 1 line if MATLAB was cooler.

h = (J>=I) .* R  + (J<I).* C;

end

