function Ht = transpose(H)
% PLAIN (non-conjugate) TRANSPOSE FOR HSS MATRIX -- the ".'" operator.
% NOTE: previously conjugated H.D here (via ') while odt() below only
% ever plain-transposed Z/Y/lrcomponent (via .') -- an inconsistent mix
% that happened to cancel out for real-valued matrices (where ' and .'
% agree) but produced a wrong result for complex ones. See ctranspose.m
% for the conjugating counterpart bound to the "'" operator.
if H.isleaf
    Ht = H;
    Ht.sz = [H.sz(2) H.sz(1)];
    Ht.Ir = H.Ic;
    Ht.Ic = H.Ir;
    Ht.D = (H.D).';
    return
end
Ht = H;
Ht.A11 = transpose(H.A11);
Ht.A12 = odt(H.A21);
Ht.A21 = odt(H.A12);
Ht.A22 = transpose(H.A22);
Ht.sz = [H.sz(2) H.sz(1)];
Ht.Ir = H.Ic;
Ht.Ic = H.Ir;
end



function Ht = odt(H)
% off diagonal transpose
Ht = H;

Ht.sz = [H.sz(2) H.sz(1)];
Ht.Ir = H.Ic;
Ht.Ic = H.Ir;
Ht.rowtreeindex = H.coltreeindex;
Ht.coltreeindex = H.rowtreeindex;
Ht.lowrankrows = H.lowrankcols;
Ht.lowrankcols = H.lowrankrows;
Ht.Z = (H.Y).';
Ht.lrcomponent = (H.lrcomponent).';
Ht.Y = (H.Z).';
end