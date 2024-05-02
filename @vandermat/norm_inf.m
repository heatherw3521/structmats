function h = norm_inf(V)
%NORM_INF Efficiently computes the infinity-norm for Vandemonde matrices. 
%   The inf-norm is the maximum absolute row sum.

    % Get dimensions of V
    m = size(V,1); % number of rows
    n = size(V,2); % number of columns
    x = V.x; % Second column of V, this determines the geometric progressions

    % Use geometric sum formula
    geo = @(r, n) (1-r.^n)./(1-r); % n-term sum, starting with 1

    % Formula fails for r=1, deal with this
    ax = abs(x);
    r = ax(ax~=1);
    rowsums = geo(r, n);

    h = max(n, max(rowsums));
end