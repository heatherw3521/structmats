function Ht = transpose(H)
% CONJUGATE TRANSPOSE FOR HSS MATRIX
if H.isleaf
    Ht = H;
    Ht.size = [H.size(2) H.size(1)];
    Ht.Ir = H.Ic;
    Ht.Ic = H.Ir;
    Ht.D = (H.D)';
    return
end
Ht = H;
Ht.A11 = transpose(H.A11);
Ht.A12 = odt(H.A21);
Ht.A21 = odt(H.A12);
Ht.A22 = transpose(H.A22);
Ht.size = [H.size(2) H.size(1)];
Ht.Ir = H.Ic;
Ht.Ic = H.Ir;
end



function Ht = odt(H)
% off diagonal transpose
Ht = H;

Ht.size = [H.size(2) H.size(1)];
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