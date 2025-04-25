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
    % a11rows = H.A11.Ir(1):H.A11.Ir(2);
    % a11cols = H.A11.Ic(1):H.A11.Ic(2);
    % a12rows = H.A12.Ir(1):H.A12.Ir(2);
    % a12cols = H.A12.Ic(1):H.A12.Ic(2);
    % a21rows = H.A21.Ir(1):H.A21.Ir(2);
    % a21cols = H.A21.Ic(1):H.A21.Ic(2);
    % a22rows = H.A22.Ir(1):H.A22.Ir(2);
    % a22cols = H.A22.Ic(1):H.A22.Ic(2);
    rows1 = H.A11.Ir(1):H.A11.Ir(2);
    rows2 = H.A22.Ir(1):H.A22.Ir(2);
    cols1 = H.A11.Ic(1):H.A11.Ic(2);
    cols2 = H.A22.Ic(1):H.A22.Ic(2);
    
    intr1 = intersect(xind,rows1);
    intr2 = intersect(xind,rows2);
    intc1 = intersect(yind,cols1);
    intc2 = intersect(yind,cols2);

    if isempty(intr1)==0 && isempty(intc1)==0
        s11 = extract(H.A11, intersect(xind,rows1), intersect(yind, cols1), xshift, yshift);
    else
        s11 = [];
    end

    if isempty(intr1)==0 && isempty(intc2)==0
        s12f = blockbuilder(H.A12,H);
        s12 = s12f(intersect(xind,rows1)-xshift,intersect(yind,cols2)-H.A11.Ic(2));
    else
        s12 = [];
    end

    if isempty(intr2)==0 && isempty(intc1)==0
        s21f = blockbuilder(H.A21,H);
        s21 = s21f(intersect(xind,rows2)-H.A11.Ir(2),intersect(yind,cols1)-yshift);
    else
        s21 = [];
    end

    if isempty(intr2)==0 && isempty(intc2)==0
        s22 = extract(H.A22, intersect(xind,rows2), intersect(yind, cols2), H.A11.Ir(2), H.A11.Ic(2));
    else
        s22 = [];
    end

    % s11 = extract(H.A11,intersect(xind,a11rows-H.Ir(1)+1),intersect(yind,a11cols-H.Ic(1)+1));
    % s12f = blockbuilder(H.A12,H);
    % s12 = s12f(intersect(xind,a12rows-H.Ir(1)+1),intersect(yind,a12cols-H.A11.Ic(2)));
    % s21f = blockbuilder(H.A21,H);
    % s21 = s21f(intersect(xind,a21rows-H.A11.Ir(2)),intersect(yind,a21cols-H.Ic(1)+1));
    % s22 = extract(H.A22,intersect(xind,a22rows-H.A11.Ir(2)),intersect(yind,a22cols-H.A11.Ic(2)));
    s = [s11, s12; s21, s22];

end

end
