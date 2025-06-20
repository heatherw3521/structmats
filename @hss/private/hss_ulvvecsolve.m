function x = hss_ulvvecsolve(H,b)
% solves Hx = b for x with H an HSS matrix

% at the leaf we dense solve
if H.isleaf
    x = H.D\b;
    return
end

HULV = H;

[HULV,QsT,SsT,bnew,wbar,rows,cols] = down_recurse(HULV,b,zeros(H.size(2),1),{},{},{},{});

bnewnew = bnew-HULV*wbar;
bnewnewnew = bnewnew(~wbar);
HULV = discard(HULV);
HULVcopy = HULV;
HULVcopy = levelup(HULVcopy,1,1);
x1 = hss_ulvvecsolve(HULVcopy,bnewnewnew);
wbar(~wbar) = x1;
x = blkdiag(SsT{:})*wbar;


end

function [H, QsT, SsT, b, wbar, drows, srows] = down_recurse(H, b, wbar, QsT, SsT, drows, srows)
if H.A11.isleaf
    % compute requisite factors on the leaf level for the two off diag
    % leaves
    [H,QsT,SsT, drows,srows] = leafcompression(H, QsT, SsT, drows, srows, 1);
    [H,QsT,SsT, drows,srows] = leafcompression(H, QsT, SsT,drows, srows, 2);
    % use Q to update rhs (being careful in the order in which Q is stored)
    b(1:H.A11.size(1)) = bupdate(cell2mat(QsT(end-1)),b(1:H.A11.size(1)));
    b(H.A11.size(1)+1:end) = bupdate(cell2mat(QsT(end)),b(H.A11.size(1)+1:end));
    % use piece of dense block that has been diagonalized to direct solve
    wbar(1:H.A11.size(2)) = bsolve(H.A11.D, b(1:H.A11.size(1)), size(H.A12.Z, 2));
    wbar(H.A11.size(2)+1:end) = bsolve(H.A22.D, b(H.A11.size(1)+1:end), size(H.A21.Z,2));
else
    [H.A11,QsT,SsT, b(1:H.A11.size(1)), wbar(1:H.A11.size(2)), drows,srows] = down_recurse(H.A11, b(1:H.A11.size(1)), wbar(1:H.A11.size(2)), QsT, SsT, drows, srows);
    [H.A22,QsT,SsT, b(H.A11.size(1)+1:end), wbar(H.A11.size(2)+1:end), drows,srows] = down_recurse(H.A22, b(H.A11.size(1)+1:end), wbar(H.A11.size(2)+1:end), QsT, SsT, drows, srows);
end
end

function [Hparent, Qs, SsT, drows, srows] = leafcompression(Hparent, Qs, SsT, drows, srows, whichrow)
if whichrow == 1
    H = Hparent.A12;
elseif whichrow == 2
    H = Hparent.A21;
end
% at the leaf level
k = size(H.Z,2); % low rank 
Fleft = flip(eye(size(H.Z,1))); % large anti diag of 1's for QL decomps
Fright = flip(eye(k)); % small anti diag of 1's for QL decomps
%drows{length(drows)+1} = [1:H.size(1)-k]; % partition row into length-k and k sections
%srows{length(srows)+1} = [H.size(1)-k+1:H.size(1)];

% ZF = QR  
[Q,R] = qr(H.Z*Fright);
% Z = QRF = (QF)(FRF) -> (QF)'Z = FRF = newZ
H.Z = Fleft*R*Fright;
Q = Q*Fleft;
Qs = [Qs; Q'];

% upper right or lower left leaf
if mod(H.rowtreeindex,2) == 1
    % update parent
    Hparent.A12 = H;
    Hparent.A11.D = Q' * Hparent.A11.D;

    D = Hparent.A11.D;
    l = size(D,1);
    % make dense block lower triangular
    % D = [D1;D2]
    % want D*Z = [L'; stuff]
    % D1' = S*L -> D1*S = L'
    [S,L] = qr(D(1:l-k,:)');
    D(1:l-k,:) = L';
    % D2 = D2*S
    D(l-k+1:end,:) = D(l-k+1:end,:)*S;

    % update dense
    Hparent.A11.D = D;

    % update column
    Hparent.A21.Y = Hparent.A21.Y * S;
    SsT = [SsT; S'];
else
    % update parent
    Hparent.A21 = H;
    Hparent.A22.D = Q' * Hparent.A22.D;

    D = Hparent.A22.D;
    l = size(D,1);
    % make dense block lower triangular
    % D = [D1;D2]
    % want D*Z = [L'; stuff]
    % D1' = S*L -> D1*S = L'
    [S,L] = qr(D(1:l-k,:)');
    D(1:l-k,:) = L';
    % D2 = D2*S
    D(l-k+1:end,:) = D(l-k+1:end,:)*S;

    % update
    Hparent.A22.D = D;

    % update column
    Hparent.A12.Y = Hparent.A12.Y * S;
    %SsT = [SsT(1:end-1); S'; SsT(end)];
    SsT = [SsT; S'];
end
end

function H = discard(H)
if H.A11.isleaf
    m1 = H.A11.size(1);
    n1 = H.A11.size(2);
    m2 = H.A22.size(1);
    n2 = H.A22.size(2);
    k1 = size(H.A12.Z,2);
    k2 = size(H.A21.Z,2);
    % discard first l1-k rows
    % from off diag
    H.A12.Z = H.A12.Z(m1-k1+1:end,:);
    H.A12.Y = H.A12.Y(:,n2-k1+1:end);
    H.A12.size = [size(H.A12.Z,1), size(H.A12.Y,2)];
    % discard first l2-k columns
    % from off diag
    H.A21.Z = H.A21.Z(m2-k2+1:end,:);
    H.A21.Y = H.A21.Y(:,n1-k2+1:end);
    H.A21.size = [size(H.A21.Z,1), size(H.A21.Y,2)];
    
    
    % from diag
    bigm1 = max(m1,n1);
    bigm2 = max(m2,n2);
    H.A11.D = H.A11.D(bigm1-k1+1:end,bigm1-k2+1:end);
    H.A11.size = size(H.A11.D);
    H.A22.D = H.A22.D(bigm2-k2+1:end,bigm2-k1+1:end);
    H.A22.size = size(H.A22.D);
else
    H.A11 = discard(H.A11);
    %H.A11.size = size(H.A11.A11) + size(H.A11.A22);
    H.A22 = discard(H.A22);
    %H.A22.size = size(H.A22.A11) + size(H.A22.A22);
end
end

function [H,rstart,cstart] = levelup(H,rstart,cstart)
if H.A11.isleaf
    H.D = [H.A11.D, H.A12.Z*H.A12.lrcomponent*H.A12.Y; H.A21.Z*H.A21.lrcomponent*H.A21.Y, H.A22.D];
    H.A11 = [];
    H.A22 = [];
    H.A12 = [];
    H.A21 = [];
    H.isleaf = 1;
elseif H.A11.A11.isleaf
    [H,rstart,cstart] = onerebuild(H.A12,H,rstart,cstart);
    [H,rstart,cstart] = onerebuild(H.A21,H,rstart,cstart);
    [H,rstart,cstart] = onerebuild(H.A11,H,rstart,cstart);
    [H,rstart,cstart] = onerebuild(H.A22,H,rstart,cstart);
    H.A12.Ir = H.A11.Ir;
    H.A12.Ic = H.A22.Ic;
    H.A21.Ir = H.A22.Ir;
    H.A21.Ic = H.A11.Ic;
    H.size = H.A11.size + H.A22.size;
    H.Ir = [H.A11.Ir(1) H.A22.Ir(2)];
    H.Ic = [H.A11.Ic(1) H.A22.Ic(2)];
else
    [H.A11,rstart,cstart] = levelup(H.A11,rstart,cstart);
    H.A22 = levelup(H.A22,rstart,cstart);
end
end

function [Hparent,rstart,cstart] = onerebuild(H,Hparent,rstart,cstart)
if H.isdiag
    % create dense block on level up on the diag
    H.D = [H.A11.D, H.A12.Z*H.A12.lrcomponent*H.A12.Y; H.A21.Z*H.A21.lrcomponent*H.A21.Y,H.A22.D];
    H.A11 = [];
    H.A12 = [];
    H.A21 = [];
    H.A22 = [];
    H.size = size(H.D);
    H.isleaf = 1;
    H.Ir = [rstart, rstart + H.size(1)-1];
    H.Ic = [cstart, cstart + H.size(2)-1];
    rstart = rstart + H.size(1);
    cstart = cstart + H.size(2);
    if mod(H.rowtreeindex,2) ==1
        Hparent.A11 = H;
    else
        Hparent.A22 = H;
    end
elseif mod(H.rowtreeindex,2) ==1
    % offdiag layer to be a leaf layer
    Z1 = Hparent.A11.A12.Z;
    Z2 = Hparent.A11.A21.Z;
    Hparent.A12.Z = blkdiag(Z1,Z2)*Hparent.A12.Z;
    
    Y1 = Hparent.A22.A21.Y;
    Y2 = Hparent.A22.A12.Y;
    Hparent.A12.Y = Hparent.A12.Y*blkdiag(Y1,Y2);

    Hparent.A12.size = [size(Hparent.A12.Z,1),size(Hparent.A12.Y,2)];

elseif mod(H.rowtreeindex,2) == 0
    % offdiag layer to be leaf layer
    Z1 = Hparent.A22.A12.Z;
    Z2 = Hparent.A22.A21.Z;
    Hparent.A21.Z = blkdiag(Z1,Z2)*Hparent.A21.Z;
    
    Y1 = Hparent.A11.A21.Y;
    Y2 = Hparent.A11.A12.Y;
    Hparent.A21.Y = Hparent.A21.Y*blkdiag(Y1,Y2);

    Hparent.A21.size = [size(Hparent.A21.Z,1),size(Hparent.A21.Y,2)];
end

end

function b = bupdate(Q, b)
b = Q*b;
end

function wbar = bsolve(D, b, k)
[m,n] = size(D);
bigm = max([m,n]);
w = D(1:bigm-k,1:bigm-k)\b(1:bigm-k);
wbar = zeros(n,1);
wbar(1:size(w,1)) = w;
end