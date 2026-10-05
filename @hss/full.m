function D = full(H)
%FULL  Dense matrix represented by an HSS matrix (O(m n) memory).
if H.isleaf
    D = leafbuild(H);
    return
end
if H.isdiag
    D = [blockbuilder(H.A11,H),blockbuilder(H.A12,H);blockbuilder(H.A21,H),blockbuilder(H.A22,H)];
else
    error('hss:full:node', 'full needs a diagonal node of an HSS tree.')
end
end

