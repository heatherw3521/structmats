function s = extract(H, xind, yind, xshift, yshift)
%EXTRACT 
%xind = cell2mat(xind)
%yind = cell2mat(yind)
if isempty(xind)
    s = [];
elseif isempty(yind)
    s = [];
elseif H.isleaf
    if H.isdiag
        sf = H.D;
        %s = sf(xind-H.Ir(1)+1,yind-H.Ic(1)+1);
    else
        sf = H.Z*H.lrcomponent*H.Y;
        %s = sf(xind-H.Ir(1)+1,yind-H.Ic(1)+1);
    end
    s = sf(xind-xshift,yind-yshift);
else
    rows1 = H.A11.Ir(1):H.A11.Ir(2);
    rows2 = H.A22.Ir(1):H.A22.Ir(2);
    cols1 = H.A11.Ic(1):H.A11.Ic(2);
    cols2 = H.A22.Ic(1):H.A22.Ic(2);

    % NOTE: previously split xind/yind across quadrants with intersect(),
    % which sorts its result ascending and drops duplicates -- unlike
    % A(xind,yind), which preserves the caller's order and repeats. So
    % H([5 2 8],:) silently came back sorted/deduped instead of matching
    % A([5 2 8],:). ismember()+logical indexing below preserves xind's/
    % yind's original order and duplicates through the split, and each
    % quadrant's result is scattered back into its original positions in
    % s (rather than block-concatenated), so interleaved requests that mix
    % rows1/rows2 (or cols1/cols2) in arbitrary order still land correctly.
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
