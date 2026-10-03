function Ht = ctranspose(H)
% CONJUGATE TRANSPOSE FOR HSS MATRIX -- the "'" operator.
%
% Previously unimplemented: with no ctranspose.m in this folder, H' on a
% scalar (1x1) hss object hit MATLAB's default array-level ctranspose,
% which is a no-op for a 1x1 object -- so H' silently returned H itself,
% completely unchanged, rather than erroring or transposing. This mirrors
% transpose.m (the ".'" operator) structurally, but conjugates every
% stored numeric factor (D, Z, Y, lrcomponent) instead of merely
% rearranging them, matching real MATLAB ' semantics.
if H.isleaf
    Ht = H;
    Ht.size = [H.size(2) H.size(1)];
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
Ht.size = [H.size(2) H.size(1)];
Ht.Ir = H.Ic;
Ht.Ic = H.Ir;
end



function Ht = odt(H)
% off diagonal conjugate transpose
Ht = H;

Ht.size = [H.size(2) H.size(1)];
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
