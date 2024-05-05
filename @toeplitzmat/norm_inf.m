function h = norm_inf(T)
%NORM_INF Efficiently computes the infinity-norm for Toeplitz matrices. 
%   The inf-norm is the maximum absolute row sum.

    % Get dimensions of T
    m = size(T,1);
    n = size(T,2);

    rc = [T.tr(n:-1:1).';T.tc(2:m-n+1)];
    rc = rc(1:m-1); % Rightmost column of T, excluding last element
    lc = [T.tc(2:end)]; % Leftmost column of T, excluding first element
    d = [0 ; cumsum(abs(lc)-abs(rc))]; % Difference of consecutive row sums
    rowsums = sum(abs(T.tr)) + d;
    h = max(rowsums);
end