function h = norm_inf(T)
%NORM_INF Efficiently computes the infinity-norm for Toeplitz matrices. 
%   The inf-norm is the maximum absolute row sum.

    % Get dimensions of T
    m = size(T,1);
    n = size(T,2);

    rs = [T.tr(n:-1:1).';T.tc(2:m-n+1)];
    rs = rs(1:m-1); % Rightmost column of T, excluding last element
    ls = [T.tc(2:end)]; % Leftmost column of T, excluding first element
    d = [0 ; cumsum(abs(ls)-abs(rs))]; % Difference of consecutive row sums
    rowsums = sum(abs(T.tr)) + d;
    h = max(rowsums);
end

