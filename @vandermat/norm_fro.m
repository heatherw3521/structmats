function h = norm_fro(V)
%NORM_FRO Efficiently computes the frobenius norm of Vandermonde matrices.
%   The Frobenius norm is defined as the square-root of the sum of 
%   the squared entries.

    m = size(V,1); % number of rows
    n = size(V,2); % number of columns
    x = V.x; % Second column of V, this determines the geometric progressions
        
    % Use geometric sum formula
    geo = @(r, n) (1-r.^n)./(1-r); % n-term sum, starting with 1

    % Formula fails for r=1, deal with this
    i1 = (abs(x)==1);
    num1 = nnz(i1);

    % Compute norm
    r = abs(x(~i1)).^2;
    h = sqrt(sum(geo(r, n)) + n*num1);
end

