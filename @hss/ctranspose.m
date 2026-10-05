function Ht = ctranspose(H)
%CTRANSPOSE  Conjugate transpose H' of an HSS matrix.
%   Swaps the row and column generators and conjugate-transposes every
%   stored block (D, Z, Y, lrcomponent); the result is HSS on the
%   transposed tree. See also transpose.
if H.isleaf
    Ht = H;
    Ht.sz = [H.sz(2) H.sz(1)];
    Ht.Ir = H.Ic;
    Ht.Ic = H.Ir;
    Ht.D = (H.D)';
    return
end
Ht = H;
Ht.A11 = ctranspose(H.A11);
Ht.A12 = odt(H.A21);
Ht.A21 = odt(H.A12);
Ht.A22 = ctranspose(H.A22);
Ht.sz = [H.sz(2) H.sz(1)];
Ht.Ir = H.Ic;
Ht.Ic = H.Ir;
end



function Ht = odt(H)
% off diagonal conjugate transpose
Ht = H;

Ht.sz = [H.sz(2) H.sz(1)];
Ht.Ir = H.Ic;
Ht.Ic = H.Ir;
Ht.rowtreeindex = H.coltreeindex;
Ht.coltreeindex = H.rowtreeindex;
Ht.lowrankrows = H.lowrankcols;
Ht.lowrankcols = H.lowrankrows;
Ht.Z = (H.Y)';
Ht.lrcomponent = (H.lrcomponent)';
Ht.Y = (H.Z)';
end
