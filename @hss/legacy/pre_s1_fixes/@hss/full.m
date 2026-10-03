function D = full(H)
% full HSS
% NOTE: previously had no base case for H.isleaf -- full() on a
% single-leaf HSS matrix (whole matrix <= blocksize) crashed trying to
% recurse into H.A11 etc., which don't exist for a leaf.
if H.isleaf
    D = leafbuild(H);
    return
end
if H.isdiag
    D = [blockbuilder(H.A11,H),blockbuilder(H.A12,H);blockbuilder(H.A21,H),blockbuilder(H.A22,H)];
else
    error('i messed this up, need to fix')
end
end

