function F = hss_cachedfactor(H)
%HSS_CACHEDFACTOR  ULV factors of H (square or wide), from H's cache if there.
%   The first call factors H and stores the result in H.factorcache (a handle,
%   so the caller's H keeps it); later calls return it at no cost.
C = H.factorcache;
if ~isempty(C) && ~isempty(C.ulv)
    F = C.ulv;
    return
end
F = hss_ulvminnormsolve(H);
if ~isempty(C)
    C.ulv = F;
end
end
