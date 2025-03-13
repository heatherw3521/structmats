function h = norm_inf(V)
%NORM_INF Efficiently computes the infinity-norm for Vandemonde matrices. 
%   The inf-norm is the maximum absolute row sum.

    n = size(V,2); % number of columns

    r = max(abs(V.x)); % Maximal geometric ratio
    if r == 1
        h = n;
    else
        h = (1-r.^n)./(1-r); % Geometric sum
    end
end