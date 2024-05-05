function h = abs(T)
%ABS Elementwise absolute value for Diagonally Structured matrices

% little hack
% we can do logic on numerical data too
ftest = @(x) true;
h = diag_apply(@(x,y) abs(x), ftest, "| |", T, 1);
end

