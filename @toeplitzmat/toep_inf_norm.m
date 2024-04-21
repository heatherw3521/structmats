function h = toep_inf_norm(T)
%TOEP_INF_NORM Quickly computes the infinity-norm for Toeplitz matrices. 
%   The inf-norm is the maximum absolute row sum.

    rs = [flip(cumsum(T.tr(2:end))) 0]; % super-diag partial sums
    cs = [0 cumsum(T.tc(2:end))']; % sub-dia partial sums
    h = T.tc(1) + max(rs+cs);
end

