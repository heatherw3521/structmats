function H = hss_scale(H, s)
%HSS_SCALE  s*H for a scalar s: scales every dense diagonal block D and every
%   coupling matrix (lrcomponent); the bases are unchanged.
if H.isdiag
    if H.isleaf
        H.D = s*H.D;
    else
        H.A11 = hss_scale(H.A11, s);
        H.A22 = hss_scale(H.A22, s);
        H.A12.lrcomponent = s*H.A12.lrcomponent;
        H.A21.lrcomponent = s*H.A21.lrcomponent;
    end
end
end
