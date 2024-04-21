function h = norm_inf(T)
%NORM_INF Quickly computes the infinity-norm for Toeplitz matrices. 
%   The inf-norm is the maximum absolute row sum.

    rs = [flip(cumsum(abs(T.tr(2:end)))) 0]; % super-diag partial sums
    cs = [0 cumsum(abs(T.tc(2:end)))']; % sub-diag partial sums
    h = T.tc(1) + max(rs+cs);
end

