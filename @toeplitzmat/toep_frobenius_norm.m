function h = toep_frobenius_norm(T)
%TOEP_FROBENIUS_NORM Efficiently computes the frobenius norm of a Toeplitz
%   matrix
    m = size(T,1);
    n = size(T,2);
    subsum = sum((m:-1:1)' .* (T.tc).^2); % Contribution on & below the diag
    supsum = sum((n-1:-1:1) .* (T.tr(2:end)).^2); % Contribution on & above diag
    h = sqrt(subsum + supsum);
end

