function s = extract(H, xind, yind, xshift, yshift)
%EXTRACT  S = H(xind, yind) for global index vectors (any order, repeats
%   allowed). xshift/yshift are the row/column offsets of H in the full
%   matrix (0 at the root). Called by subsref after the indices are checked.
if isempty(xind)
    s = [];
elseif isempty(yind)
    s = [];
elseif H.isleaf
    if H.isdiag
        sf = H.D;
    else
        sf = H.Z*H.lrcomponent*H.Y;
    end
    s = sf(xind-xshift,yind-yshift);
else
    rows1 = H.A11.Ir(1):H.A11.Ir(2);
    rows2 = H.A22.Ir(1):H.A22.Ir(2);
    cols1 = H.A11.Ic(1):H.A11.Ic(2);
    cols2 = H.A22.Ic(1):H.A22.Ic(2);

    % split the requested rows/columns between the two halves with masks,
    % so that the caller's order and repeated indices are kept, and scatter
    % each quadrant back into its positions in s
    rmask1 = ismember(xind, rows1);
    rmask2 = ismember(xind, rows2);
    cmask1 = ismember(yind, cols1);
    cmask2 = ismember(yind, cols2);

    s = zeros(numel(xind), numel(yind));

    if any(rmask1) && any(cmask1)
        s(rmask1, cmask1) = extract(H.A11, xind(rmask1), yind(cmask1), xshift, yshift);
    end

    if any(rmask1) && any(cmask2)
        s12f = blockbuilder(H.A12,H);
        s(rmask1, cmask2) = s12f(xind(rmask1)-xshift, yind(cmask2)-H.A11.Ic(2));
    end

    if any(rmask2) && any(cmask1)
        s21f = blockbuilder(H.A21,H);
        s(rmask2, cmask1) = s21f(xind(rmask2)-H.A11.Ir(2), yind(cmask1)-yshift);
    end

    if any(rmask2) && any(cmask2)
        s(rmask2, cmask2) = extract(H.A22, xind(rmask2), yind(cmask2), H.A11.Ir(2), H.A11.Ic(2));
    end

end

end
