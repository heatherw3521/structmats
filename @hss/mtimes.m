function H = mtimes(H1,H2)
%MTIMES  Products with an HSS matrix H.
%   H*X, X*H     dense X (one or more columns / rows)
%   s*H, H*s     scalar s (the result is HSS)
%   H1*H2        two HSS matrices on conforming trees; the result is HSS,
%                recompressed with compress(., 1e-12) since the ranks add
if isa(H1,'hss'), hss_assertroot(H1, 'H*X'); end
if isa(H2,'hss'), hss_assertroot(H2, 'X*H'); end

if isa(H1,'hss')
    if isa(H2,'hss')
        H = compress(hss_matmat(H1,H2), 1e-12);
    elseif isscalar(H2)
        H = hss_scale(H1, H2);
    else
        if size(H2, 1) ~= H1.sz(2)
            error('hss:mtimes:dimension', 'Inner matrix dimensions must agree: X has %d rows; H has %d columns.', size(H2, 1), H1.sz(2));
        end
        H = hss_matvec(H1,H2);
    end
else
    if isscalar(H1)
        H = hss_scale(H2, H1);
    else
        % X*H = (H.' * X.').'
        if size(H1, 2) ~= H2.sz(1)
            error('hss:mtimes:dimension', 'Inner matrix dimensions must agree: X has %d columns; H has %d rows.', size(H1, 2), H2.sz(1));
        end
        H = hss_matvec(transpose(H2), H1.').';
    end
end

end
