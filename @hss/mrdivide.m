function X = mrdivide(B, H)
%MRDIVIDE  X = B/H, the solution of X*H = B, for a dense B and an HSS H.
%   square H : X = B*inv(H)
%   wide H   : least squares, argmin ||X*H - B|| (H of full row rank)
%   tall H   : the minimum-norm solution of X*H = B
%   B/H reuses the ULV factors H keeps for H\b (square, wide) or the
%   factors of H' (tall). H/s for a scalar s is H*(1/s).
if isa(B, 'hss')
    if isscalar(H) && ~isa(H, 'hss')
        X = hss_scale(B, 1/H);
        return
    end
    error('hss:mrdivide:operand', 'H/X is only supported for a scalar X; for X/H the HSS matrix must be on the right.');
end
hss_assertroot(H, 'B/H');
if size(B, 2) ~= H.sz(2)
    error('hss:mrdivide:dimension', 'B has %d columns; H has %d columns.', size(B, 2), H.sz(2));
end
if H.sz(1) > H.sz(2)
    % X*H = B underdetermined: minimum-norm X' = pinv(H')*B'
    X = hss_ulvminnormsolve(hss_cachedadjoint(H), B')';
else
    % X' = pinv(H)'*B', the adjoint of the solve with H's own factors
    F = hss_cachedfactor(H);
    if ~hss_ulv_fullrank(F)
        error('hss:mrdivide:rankDeficient', ['H does not have full row rank, so B/H ' ...
            'is not unique; use tikhonov(H'', B'', lambda)'' instead.'])
    end
    X = hss_ulv_adjoint(F, B')';
end
end
