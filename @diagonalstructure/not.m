function h = not(T)
%~  Pointwise logical not for DIAGONALSTRUCTURE objects
%   T is a DIAGONALSTRUCTURE
%   behaves the same as ~ for matrices in MATLAB

% little hack because I'm too lazy to change toepcompare
% we can do logic on numerical data too
ftest = @(x) true;
h = diag_apply(@(x,y) ~x, ftest, "~", T, 1);
end