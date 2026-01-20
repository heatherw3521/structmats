function [x,H] = hss_udsolve(H,b,options)
% solves Hx = b for x with H an HSS matrix and m<n
arguments
    H
    b
    options.minnormflag = 1;
    options.storageflag = 1;
end

if H.size(1)>H.size(2)
    error("Underdetermined Solve (UD) assumes more columns than rows!")
end

if H.size(1) == H.size(2)
    x = hss_ulvvecsolve(H,b);
    return
end
if H.isleaf
    x = lsqminnorm(H.D,b);
    return
end

b_init = b;
H_ud = H;

x_init = zeros(H.size(2),1);
Ss = {};

[H_ud,b,x,Ss,r_inds,H] = down_recurse(H_ud,b,x_init,Ss,H);

bnew = b - H_ud*x;

%binds = ~ismembertol(b,H_ud*x);
binds = ~ismembertol(bnew,zeros(size(b)));
bnew = bnew(binds);
%xnew = x(~x);

H_udnew = discard(H_ud);
H_udnew = levelup(H_udnew,1,1);

x1 = hss_udsolve(H_udnew,bnew);
%x(~x) = x1;
x(r_inds) = x1;
x = (blkdiag(Ss{:}))*x;

if options.minnormflag
    [x,H] = minxnorm(H,x,b_init,options.storageflag);
end

end


function [H,b,x,Ss,r_inds,H_orig] = down_recurse(H,b,x,Ss,H_orig)

% the leaf level we do the following:
% eliminate rows on the off-diags via QR
% apply that transform to the diag and to the rhs
% antidiagonalize dense component via QR
% apply that transform to the columns
% store those transforms to reapply to the solution vector

if H.A11.isleaf
    % find the transforms to eliminate the first n-k rows of each leaf
    % if it already exists for the given H apply it instead
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

    % apply inverse of transforms to the rhs 
    b(1:H.A11.size(1)) = Q1'*b(1:H.A11.size(1));
    b(H.A11.size(1)+1:end) = Q2'*b(H.A11.size(1)+1:end);
    % apply transforms to the dense component
    H.A11.D = Q1'*H.A11.D;
    H.A22.D = Q2'*H.A22.D;

    % via QR find transform to "triangularize" dense component (0 above shifted
    % diagonal to opposite corner)
    % again, if it already exists for H apply it instead
    
    if H_orig.A11.S
        S1 = H_orig.A11.S;
        H.A11.D = H.A11.D*S1;
    else
        [S1,H.A11] = diagonalcompression(H.A11);
        H_orig.A11.S = S1;
    end
    
    if H_orig.A11.S
        S2 = H_orig.A22.S;
        H.A22.D = H.A22.D*S2;
    else
        [S2,H.A22] = diagonalcompression(H.A22);
        H_orig.A22.S = S2;
    end

    % apply transform to column space
    H.A21.Y = H.A21.Y*S1;
    H.A12.Y = H.A12.Y*S2;
    % store S to apply S' to solution x at the end
    Ss = [Ss;S1;S2];

    % solve x components in overdetermined system with (l1-k)x(l2-k) block
    % in diag blocks

    % A11
    l1_1 = H.A11.size(1);
    l2_1 = H.A11.size(2);
    k_1 = size(H.A12.Z,2);
    if l1_1-k_1 == 0
        error("Off-diagonal components are not compressable. Solve directly with the full matrix.")
    end
    x(1:(l2_1-k_1)) = lsqminnorm(H.A11.D(1:(l1_1-k_1),1:(l2_1-k_1)), b(1:(l1_1-k_1)));

    % A22
    l1_2 = H.A22.size(1);
    l2_2 = H.A22.size(2);
    k_2 = size(H.A21.Z,2);
    if l1_2-k_2 == 0
        error("Off-diagonal components are not compressable. Solve directly with the full matrix.")
    end

    x((l2_1+1):(end-k_2)) = lsqminnorm(H.A22.D(1:(l1_2-k_2),1:(l2_2-k_2)), b((l1_1+1):(end-k_2)));
    r_inds = [l2_1-k_1+1:l2_1,size(x,1)-k_2+1:size(x,1)];

% on the non-leaf level we recurse down upper left and lower right    
else
    [H.A11,b(1:H.A11.size(1)),x(1:H.A11.size(2)),Ss,r_inds1,H_orig.A11] = down_recurse(H.A11,b(1:H.A11.size(1)),x(1:H.A11.size(2)),Ss,H_orig.A11);
    [H.A22,b(H.A11.size(1)+1:end),x(H.A11.size(2)+1:end),Ss,r_inds2,H_orig.A22] = down_recurse(H.A22,b(H.A11.size(1)+1:end),x(H.A11.size(2)+1:end),Ss,H_orig.A22);
    r_inds2 = r_inds2+H.A11.size(2);
    r_inds = [r_inds1,r_inds2];
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

%Fleft = flip(eye(l1));
%Fright = flip(eye(l2));

D1 = H.D(:,1:l2-l1);
D2 = H.D(:,l2-l1+1:end);

[Q,R] = qr(D2');

Q = blkdiag(eye(abs(l1-l2)),Q);

H.D = [D1 R'];

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

function [x,H] = minxnorm1(H,x,b,storageflag)
%% min norm
if H.minnormQ
    Q = H.minnormQ;
else
    diff = abs(H.size(1)-H.size(2));
    %iters = max(floor(diff/3),100);
    iters = diff;
    N = size(x,1);
    Null = zeros(N,iters);
    for i=1:iters
        v = randn(N,1);
        bk = b+H*v; %nlogn
        % xk = H\bk; %
        xk = hss_udsolve(H,bk,minnormflag = 0);
        if norm(H*xk-bk)>1e-6
            break
        end
        z = xk-x-v;
        Null(:,i) = z;
    end
    
    [Q,~] = qr(Null,'econ');
    if storageflag
        H.minnormQ = Q;
    end
end
x = x-Q*Q'*x;
end

function [x,H] = minxnorm(H,x,b,storageflag)
%% minimum-norm projection helper

if H.minnormQ
    Q = H.minnormQ;
else
    m = H.size(1);
    n = H.size(2);
    diff = abs(m - n);   % proxy for nullity
    N = n;

    % Heuristic split
    %use_range = (m < n) && (m < diff);  % r < d regime

    global use_range
    if ~use_range
        %% NULLSPACE CONSTRUCTION
        iters = diff;
        Null = zeros(N,iters);

        for i = 1:iters
            v = randn(N,1);
            bk = b + H*v;                 % A*v
            xk = hss_udsolve(H,bk,minnormflag = 0);

            if norm(H*xk - bk) > 1e-6
                Null = Null(:,1:i-1);
                break
            end

            z = xk - x - v;               % nullspace vector
            Null(:,i) = z;
        end

        [Q,~] = qr(Null,'econ');
        nr = 1;

    else
        %% RANGE(A^T) CONSTRUCTION (cheap when r << d)
        % oversample = 10;
        % iters = min(m, diff) + oversample;
        %iters = min(m,diff);

        V = randn(m,m);
        Y = H.'*V;
        % Y = zeros(n,iters);
        % 
        % for i = 1:iters
        %     Y(:,i) = H' * V(:,i);         % A^T * v
        % end

        [Q,~] = qr(Y,'econ');
        nr = 0;
    end

    if storageflag
        H.minnormQ = Q;
        H.nullorrange = nr;
    end
end

%% Final projection
if H.nullorrange == 1
    % Q spans nullspace
    x = x - Q*(Q'*x);
else
    % Q spans range(A^T)
    x = Q*(Q'*x);
end
end

