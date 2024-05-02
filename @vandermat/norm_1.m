function h = norm_1(V)
%NORM_INF Efficiently computes the 1-norm for Vandemonde matrices. 
%   The 1-norm is the maximum absolute column sum.


    % We have a convex problem on an interval! The maximum is at one of the
    % endpoints.

    % Get dimensions of V
    m = size(V,1); % number of rows
    n = size(V,2); % number of columns
    ax = abs(V.x); % Second column of V, this determines the geometric progressions

    h = max(m, sum(ax .^ (n-1)));
end