function [x,H] = hss_ulvvecsolve(H,b)
% solves Hx = b for x with H an HSS matrix

% at the leaf we dense solve
if H.isleaf
    x = H.D\b;
    return
end

HULV = H;

[HULV,Ss,bnew,wbar,H] = down_recurse(HULV,b,zeros(H.size(2),1),{},H);

bnewnew = bnew-HULV*wbar;
bnewnewnew = bnewnew(~wbar); % THIS LINE NEEDS TO BE READDRESSED BC I THINK ITS WRONG FOR RECT CASE
HULV = discard(HULV);
%Hdisc = discard(HULV);
%b_reduced = bnewnew(Hdisc.Ir);
%wbar_reduced = wbar(Hdisc.Ic);

HULVcopy = HULV;
HULVcopy = levelup(HULVcopy,1,1);
% HULVcopy = levelestup(HULVcopy,1,1);
x1 = hss_ulvvecsolve(HULVcopy,bnewnewnew);
wbar(~wbar) = x1;
x = (blkdiag(Ss{:}))*wbar;


end

function [H, Ss, b, wbar,H_orig] = down_recurse(H, b, wbar, Ss, H_orig)

if H.A11.isleaf
    % compute requisite factors on the leaf level for the two off diag
    % leaves
    if H_orig.A12.Q
        Q1 = H_orig.A12.Q;
        H.A12.Z = Q1'*H.A12.Z;
    else
        [Q1,H.A12] = leafcompression(H.A12);
        H_orig.A12.Q = Q1;
    end

    if H_orig.A21.Q
        Q2 = H_orig.A21.Q;
        H.A21.Z = Q2'*H.A21.Z;
    else
        [Q2,H.A21] = leafcompression(H.A21);
        H_orig.A21.Q = Q2;
    end

    b(1:H.A11.size(1)) = Q1'*b(1:H.A11.size(1));
    b(H.A11.size(1)+1:end) = Q2'*b(H.A11.size(1)+1:end);
    H.A11.D = Q1'*H.A11.D;
    H.A22.D = Q2'*H.A22.D;

    if H_orig.A11.S
        S1 = H_orig.A11.S;
        H.A11.D = H.A11.D*S1;
    else
        [S1,H.A11] = diagonalcompression(H.A11);
        H_orig.A11.S = S1;
    end
    
    if H_orig.A22.S
        S2 = H_orig.A22.S;
        H.A22.D = H.A22.D*S2;
    else
        [S2,H.A22] = diagonalcompression(H.A22);
        H_orig.A22.S = S2;
    end

    H.A21.Y = H.A21.Y*S1;
    H.A12.Y = H.A12.Y*S2;
    Ss = [Ss;S1;S2];
    % [H,Qs,Ss, drows,srows] = leafcompression(H, Qs, Ss, drows, srows, 1);
    % [H,Qs,Ss, drows,srows] = leafcompression(H, Qs, Ss,drows, srows, 2);
    % use Q to update rhs (being careful in the order in which Q is stored)
    % if length(b) ~= H.A11.size(1) + H.A22.size(1)
    %     disp([length(b), H.A11.size(1) + H.A22.size(1)])
    % end
    % if length(wbar) ~= H.A11.size(2) + H.A22.size(2)
    %     disp([length(wbar), H.A11.size(2) + H.A22.size(2)])
    % end
    % b(1:H.A11.size(1)) = bupdate(cell2mat(Qs(end-1)),b(1:H.A11.size(1)));
    % b(H.A11.size(1)+1:end) = bupdate(cell2mat(Qs(end)),b(H.A11.size(1)+1:end));
    % use piece of dense block that has been triangularized to direct solve


    % FOR A 2 LEVEL HSS AT THIS POINT WHEN I RUN
    %xtruetest = (full(H)*(blkdiag(SsT{end-1:end})))\b %I SHOULD GET XTRUE (?) AND I DON'T
    

    wbar(1:H.A11.size(2)) = bsolve(H.A11.D, b(1:H.A11.size(1)), size(H.A12.Z, 2));
    wbar(H.A11.size(2)+1:end) = bsolve(H.A22.D, b(H.A11.size(1)+1:end), size(H.A21.Z,2));

    %(blkdiag(SsT{:}))'*wbar
else
    [H.A11,Ss, b(1:H.A11.size(1)), wbar(1:H.A11.size(2)),H_orig.A11] = down_recurse(H.A11, b(1:H.A11.size(1)), wbar(1:H.A11.size(2)), Ss, H_orig.A11);
    [H.A22,Ss, b(H.A11.size(1)+1:end), wbar(H.A11.size(2)+1:end),H_orig.A22] = down_recurse(H.A22, b(H.A11.size(1)+1:end), wbar(H.A11.size(2)+1:end), Ss, H_orig.A22);
end
end

function [Q,H] = leafcompression(H)
% Z = QRF = (QF)(FRF) -> (QF)'Z = FRF = newZ
Fleft = flip(eye(size(H.Z,1))); % large anti diag of 1s
Fright = flip(eye(size(H.Z,2))); % small anti diag of 1s
[Q,R] = qr(H.Z*Fright);
H.Z = Fleft * (R) * Fright;
Q = (Q*Fleft);
end


function [Q,H] = diagonalcompression(H)
l1 = H.size(1);
l2 = H.size(2);
D1 = H.D(:,1:l2-l1);
D2 = H.D(:,l2-l1+1:end);
[Q,R] = qr(D2');
Q = blkdiag(eye(abs(l1-l2)),Q);
H.D = [D1 R'];
end

% function [Hparent, Qs, Ss, drows, srows] = leafcompression(Hparent, Qs, Ss, drows, srows, whichrow)
% if whichrow == 1
%     H = Hparent.A12;
% elseif whichrow == 2
%     H = Hparent.A21;
% end
% % at the leaf level
% k = size(H.Z,2); % low rank 
% Fleft = flip(eye(size(H.Z,1))); % large anti diag of 1's for QL decomps
% Fright = flip(eye(k)); % small anti diag of 1's for QL decomps
% %drows{length(drows)+1} = [1:H.size(1)-k]; % partition row into length-k and k sections
% %srows{length(srows)+1} = [H.size(1)-k+1:H.size(1)];
% 
% % ZF = QR  
% [Q,R] = qr(H.Z*Fright);
% % Z = QRF = (QF)(FRF) -> (QF)'Z = FRF = newZ
% H.Z = Fleft * (R) * Fright;
% %Q = Q*Fleft;
% %Qs = [Qs; Q'];
% Q = (Q*Fleft);
% Qs = [Qs; Q];
% 
% % upper right or lower left leaf
% if mod(H.rowtreeindex,2) == 1
%     % update parent
%     Hparent.A12 = H;
%     Hparent.A11.D = Q' * Hparent.A11.D;
% 
%     D = Hparent.A11.D;
%     %l = size(D,1);
%     % make dense block lower triangular
%     % D = [D1;D2]
%     % want D*Z = [L'; stuff]
%     % D1' = S*L -> D1*S = L'
% 
% 
%     %[S,L] = qr(D(1:l-k,:)');
%     %D(1:l-k,:) = L';
%     [S,L] = qr(D');
%     D = L';
%     S = S';
%     % D2 = D2*S
%     %D(l-k+1:end,:) = D(l-k+1:end,:)*S;
% 
%     % update dense
%     Hparent.A11.D = D;
% 
%     % update column
%     Hparent.A21.Y = Hparent.A21.Y * S';
%     Ss = [Ss; S];
% else
%     % update parent
%     Hparent.A21 = H;
%     Hparent.A22.D = Q' * Hparent.A22.D;
% 
%     D = Hparent.A22.D;
%     % l = size(D,1);
%     % make dense block lower triangular
%     % D = [D1;D2]
%     % want D*Z = [L'; stuff]
%     % D1' = S*L -> D1*S = L'
%     % [S,L] = qr(D(1:l-k,:)');
%     % D(1:l-k,:) = L';
%     % % D2 = D2*S
%     % D(l-k+1:end,:) = D(l-k+1:end,:)*S;
% 
% 
%     [S,L] = qr(D');
%     D = L';
%     S = S';
% 
%     % update
%     Hparent.A22.D = D;
% 
%     % update column
%     Hparent.A12.Y = Hparent.A12.Y * S';
%     %SsT = [SsT(1:end-1); S'; SsT(end)];
%     Ss = [Ss; S];
% end
% end

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
    H.A12.Z = H.A12.Z(end-k1+1:end,:);
    H.A12.Y = H.A12.Y(:,end-k2+1:end);
    H.A12.size = [size(H.A12.Z,1), size(H.A12.Y,2)];
    H.A12.isleaf = 1;
    % discard first l2-k columns
    % from off diag
    H.A21.Z = H.A21.Z(end-k2+1:end,:);
    H.A21.Y = H.A21.Y(:,end-k1+1:end);
    H.A21.size = [size(H.A21.Z,1), size(H.A21.Y,2)];
    H.A21.isleaf = 1;


    % from diag
    bigm1 = max(m1,n1);
    bigm2 = max(m2,n2);
    H.A11.D = H.A11.D(end-k1+1:end,end-k1+1:end);
    H.A11.size = size(H.A11.D);
    H.A22.D = H.A22.D(end-k2+1:end,end-k2+1:end);
    H.A22.size = size(H.A22.D);
else
    H.A11 = discard(H.A11);
    %H.A11.size = size(H.A11.A11) + size(H.A11.A22);
    H.A22 = discard(H.A22);
    %H.A22.size = size(H.A22.A11) + size(H.A22.A22);
end
end

function H = levelestup(H,rstart,cstart)
H.levelcount = H.levelcount-1;
if H.A11.isleaf
    H.isleaf = 1;
    H.D = full(H);
    %H.D = hss_matvec(H,eye(H.size(2)));
    H.A11 = [];
    H.A12 = [];
    H.A21 = [];
    H.A22 = [];
    H.size = size(H.D);
    H.Ir = [rstart rstart+H.size(1)-1];
    H.Ic = [cstart cstart+H.size(2)-1];

else
    Hcopy = H;
    H.A11 = levelestup(Hcopy.A11,rstart,cstart);
    % H.A12 = levelestup(Hcopy.A12);
    % H.A21 = levelestup(Hcopy.A21);
    H.A22 = levelestup(Hcopy.A22,rstart+H.A11.size(1),cstart+H.A11.size(2));
    H.size = H.A11.size + H.A22.size;
    H.Ir = [rstart rstart+H.size(1)-1];
    H.Ic = [cstart cstart+H.size(2)-1];
    H.A12.size = [H.A11.size(1) H.A22.size(2)];
    H.A12.Ir = H.A11.Ir;
    H.A12.Ic = H.A22.Ic;
    H.A21.size = [H.A22.size(1) H.A11.size(2)];
    H.A21.Ir = H.A22.Ir;
    H.A21.Ic = H.A11.Ic;
end

end

function [H,rstart,cstart] = levelup(H,rstart,cstart)
if H.A11.isleaf
    H.D = [H.A11.D, H.A12.Z*H.A12.lrcomponent*H.A12.Y; H.A21.Z*H.A21.lrcomponent*H.A21.Y, H.A22.D];
    %H.D = full(H);
    H.A11 = [];
    H.A22 = [];
    H.A12 = [];
    H.A21 = [];
    H.isleaf = 1;
    H.size = size(H.D);
    H.levelcount = H.levelcount-1;
elseif H.A11.A11.isleaf
    [H,rstart,cstart] = onerebuild(H.A12,H,rstart,cstart);
    [H,rstart,cstart] = onerebuild(H.A21,H,rstart,cstart);
    [H,rstart,cstart] = onerebuild(H.A11,H,rstart,cstart);
    [H,rstart,cstart] = onerebuild(H.A22,H,rstart,cstart);
    % [H,rstart,cstart] = rebuild(H,rstart,cstart);
    H.A12.Ir = H.A11.Ir;
    H.A12.Ic = H.A22.Ic;
    H.A21.Ir = H.A22.Ir;
    H.A21.Ic = H.A11.Ic;
    H.size = H.A11.size + H.A22.size;
    H.Ir = [H.A11.Ir(1) H.A22.Ir(2)];
    H.Ic = [H.A11.Ic(1) H.A22.Ic(2)];
    H.levelcount = H.levelcount-1;
else
    [H.A11,rstart,cstart] = levelup(H.A11,rstart,cstart);
    [H.A22,rstart,cstart] = levelup(H.A22,rstart,cstart);
    H.A12.Ir = H.A11.Ir;
    H.A12.Ic = H.A22.Ic;
    H.A21.Ir = H.A22.Ir;
    H.A21.Ic = H.A11.Ic;
    H.size = H.A11.size + H.A22.size;
    H.Ir = [H.A11.Ir(1) H.A22.Ir(2)];
    H.Ic = [H.A11.Ic(1) H.A22.Ic(2)];
    H.levelcount = H.levelcount-1;
end
end


function [Hn,rstart,cstart] = rebuild(Hparent,rstart,cstart)
    Hn = Hparent;

    Z1 = Hparent.A11.A12.Z;
    Z2 = Hparent.A11.A21.Z;
    Hparent.A12.Z = blkdiag(Z1,Z2)*Hparent.A12.Z;
    
    Y1 = Hparent.A22.A21.Y;
    Y2 = Hparent.A22.A12.Y;
    Hn.A12.Y = Hparent.A12.Y*blkdiag(Y1,Y2);

    Hn.A12.size = [size(Hparent.A12.Z,1),size(Hparent.A12.Y,2)];



        % offdiag layer to be leaf layer
    Z1 = Hparent.A22.A12.Z;
    Z2 = Hparent.A22.A21.Z;
    Hparent.A21.Z = blkdiag(Z1,Z2)*Hparent.A21.Z;
    
    Y1 = Hparent.A11.A21.Y;
    Y2 = Hparent.A11.A12.Y;
    Hparent.A21.Y = Hparent.A21.Y*blkdiag(Y1,Y2);

    Hparent.A21.size = [size(Hparent.A21.Z,1),size(Hparent.A21.Y,2)];



    Hold = Hparent.A11;
    
    Hold.D = [Hold.A11.D, Hold.A12.Z*Hold.A12.lrcomponent*Hold.A12.Y; Hold.A21.Z*Hold.A21.lrcomponent*Hold.A21.Y,Hold.A22.D];
    %H.D = full(H);
    Hold.A11 = [];
    Hold.A12 = [];
    Hold.A21 = [];
    Hold.A22 = [];
    Hold.size = size(Hold.D);
    Hold.isleaf = 1;
    Hold.Ir = [rstart, rstart + Hold.size(1)-1];
    Hold.Ic = [cstart, cstart + Hold.size(2)-1];
    rstart = rstart + Hold.size(1);
    cstart = cstart + Hold.size(2);

    Hn.A11 = Hold;

    Hold = Hparent.A22;

    Hold.D = [Hold.A11.D, Hold.A12.Z*Hold.A12.lrcomponent*Hold.A12.Y; Hold.A21.Z*Hold.A21.lrcomponent*Hold.A21.Y,Hold.A22.D];
    %H.D = full(H);
    Hold.A11 = [];
    Hold.A12 = [];
    Hold.A21 = [];
    Hold.A22 = [];
    Hold.size = size(Hold.D);
    Hold.isleaf = 1;
    Hold.Ir = [rstart, rstart + Hold.size(1)-1];
    Hold.Ic = [cstart, cstart + Hold.size(2)-1];
    rstart = rstart + Hold.size(1);
    cstart = cstart + Hold.size(2);

    Hn.A22 = Hold;
    
end

function [Hparent,rstart,cstart] = onerebuild(H,Hparent,rstart,cstart)
if H.isdiag
    % create dense block on level up on the diag
    H.D = [H.A11.D, H.A12.Z*H.A12.lrcomponent*H.A12.Y; H.A21.Z*H.A21.lrcomponent*H.A21.Y,H.A22.D];
    %H.D = full(H);
    H.A11 = [];
    H.A12 = [];
    H.A21 = [];
    H.A22 = [];
    H.size = size(H.D);
    H.levelcount = H.levelcount-1;
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
    H.levelcount = H.levelcount-1;
    Hparent.A12.size = [size(Hparent.A12.Z,1),size(Hparent.A12.Y,2)];

elseif mod(H.rowtreeindex,2) == 0
    % offdiag layer to be leaf layer
    Z1 = Hparent.A22.A12.Z;
    Z2 = Hparent.A22.A21.Z;
    Hparent.A21.Z = blkdiag(Z1,Z2)*Hparent.A21.Z;
    
    Y1 = Hparent.A11.A21.Y;
    Y2 = Hparent.A11.A12.Y;
    Hparent.A21.Y = Hparent.A21.Y*blkdiag(Y1,Y2);
    H.levelcount = H.levelcount-1;
    Hparent.A21.size = [size(Hparent.A21.Z,1),size(Hparent.A21.Y,2)];
end

end

function b = bupdate(Q, b)
b = Q'*b;
end

function wbar = bsolve(D, b, k)
[m,n] = size(D);
bigm = max([m,n]);
w = D(1:end-k,1:end-k)\b(1:end-k);
wbar = zeros(n,1);
wbar(1:n-k) = w;
end