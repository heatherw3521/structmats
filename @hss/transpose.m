function Ht = transpose(H)
%TRANSPOSE  Transpose H.' of an HSS matrix (no conjugation).
%   Swaps the row and column generators and transposes every stored block;
%   the result is HSS on the transposed tree. See also ctranspose.
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