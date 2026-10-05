function H = mtimes(H1,H2)
%MTIMES  H*X for an HSS matrix H and a dense X, or H1*H2 for two HSS
%   matrices on conforming trees (the product is HSS; ranks add).
if isa(H1,'hss'), hss_assertroot(H1, 'H*X'); end
if isa(H2,'hss'), hss_assertroot(H2, 'X*H'); end

if isa(H1,'hss')
    if isa(H2,'hss')
        H = hss_matmat(H1,H2);
    elseif isscalar(H2)
        error('hss:mtimes:scalar', 'Multiplying an HSS matrix by a scalar is not supported.')
    else
        H = hss_matvec(H1,H2);
    end
else
    if isscalar(H1)
        error('hss:mtimes:scalar', 'Multiplying an HSS matrix by a scalar is not supported.')
    else
        error('hss:mtimes:left', 'X*H with a dense X is not supported; use (H.''*X.'').''.')
    end
end

end

