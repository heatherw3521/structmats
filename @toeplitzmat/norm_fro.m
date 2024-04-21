function h = norm_fro(T)
%NORM_FRO Efficiently computes the frobenius norm of a Toeplitz matrices.
%   The Frobenius norm is defined as the square-root of the sum of 
%   the squared entries.

    m = size(T,1);
    n = size(T,2);
    subsum = sum((min(n,m:-1:1))' .* (T.tc).^2); % Contribution on & below the diag
    supsum = sum((min(m,n-1:-1:1)) .* (T.tr(2:end)).^2); % Contribution on & above diag
    h = sqrt(subsum + supsum);
end

