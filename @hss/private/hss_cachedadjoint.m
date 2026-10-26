function F = hss_cachedadjoint(H)
%HSS_CACHEDADJOINT  ULV factors of H' (for tall H), from H's cache if there.
%   Used by H\b and B/H for a tall H; stored in H.factorcache.adjoint.
C = H.factorcache;
if ~isempty(C) && ~isempty(C.adjoint)
    F = C.adjoint;
    return
end
F = hss_ulvminnormsolve(ctranspose(H));
if ~isempty(C)
    C.adjoint = F;
end
end
