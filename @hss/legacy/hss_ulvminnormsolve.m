function [x,H] = hss_ulvminnormsolve(H,b)
% solves Hx=b for x with H a (possibly wide, m<n) multi-level HSS matrix,
% in the minimum-norm sense: x = argmin ||x||_2 s.t. Hx=b.

% If some leaf pair lacks the slack needed for the Size Reduction step
% below (see slack_ok), the ULV reduction cannot proceed at that level no
% matter how it's driven -- the joint LQ used there needs the sibling's
% off-diagonal rank plus its own row count to fit within its own column
% count. Rather than erroring, uniformly collapse the deepest tree level
% into denser leaves (levelup already does exactly this, and is already
% used elsewhere in this file) and recheck; repeat until every leaf pair
% has enough slack, or the whole tree has collapsed to one dense block.
% This trades some avoidable densification (a whole level collapses even
% if only one pair on it needed it) for a solver that never breaks down.
while ~H.isleaf && ~slack_ok(H)
    H = levelup(H,1,1);
end

if H.isleaf
    if size(H.D,1) < size(H.D,2)
        x = pinv(H.D)*b;
    else
        x = H.D\b;
    end
    return
end

HULV = H;
[HULV,Omegas,Qs,bnew,wbar,wind,bind,H,padinfo] = down_recurse(HULV,b,zeros(H.size(2),1),{},{},[],[],H,zeros(0,2));

bnewnew = bnew - HULV*wbar;
bnewnewnew = bnewnew(setdiff(1:size(bnewnew,1),bind));

HULVcopy = discard(HULV);
HULVcopy = levelup(HULVcopy,1,1);

x1 = hss_ulvminnormsolve(HULVcopy,bnewnewnew);
wbar(setdiff(1:size(wbar,1),wind)) = x1;

% wbar lives in the size-reduced coordinate system at this point (its
% length is sum of each leaf's reduced width, not the original width).
% Q rotates it within that reduced space; Omega then lifts each leaf's
% chunk back to its own original width (zero-padding the discarded/
% null-space part) before stitching the pieces back together.
% NOTE: this file's lq()/ql() helpers return Q such that A = L*Q, i.e.
% the collected Qs/Omegas are the TRANSPOSE of the notes' Q_tau/Omega_tau
% (which satisfy A*Omega_tau = [L,0], A*Q_tau applied on the right).
% The notes' reconstruction x=Omega[x_tilde;0], x_tilde=Q*wbar therefore
% needs the transposed blocks here.
%
% Both Qs{i} and Omegas{i} are square, collected in the same left-to-right
% order as wbar's/padded's own segments -- so blkdiag(Qs{:})'*wbar (and
% the Omegas analogue) is exactly a block-diagonal operator applied to a
% vector, which never needs the combined block-diagonal matrix formed
% explicitly: apply each block to its own slice and stitch the results
% back together, same pattern hss_ulvvecsolve.m already uses for Ss.
% blkdiag(Qs{:}) alone would be length(wbar) x length(wbar) -- for a
% genuinely large matrix (tens of thousands of leaves) that is a dense
% array of a size no HSS method should ever need, and it exhausted
% available memory outright at M=20000 before this fix.
qx = zeros(size(wbar));
off = 0;
for i = 1:length(Qs)
    ni = size(Qs{i}, 1);
    idx = off+1:off+ni;
    qx(idx) = Qs{i}' * wbar(idx);
    off = off + ni;
end

padded = zeros(sum(padinfo(:,1)),1);
off_r = 0; off_p = 0;
for pidx = 1:size(padinfo,1)
    porig = padinfo(pidx,1);
    predw = padinfo(pidx,2);
    padded(off_p+1:off_p+predw) = qx(off_r+1:off_r+predw);
    off_r = off_r+predw;
    off_p = off_p+porig;
end

x = zeros(size(padded));
off = 0;
for i = 1:length(Omegas)
    ni = size(Omegas{i}, 1);
    idx = off+1:off+ni;
    x(idx) = Omegas{i}' * padded(idx);
    off = off + ni;
end

end


function [H, Omegas, Qs, b, wbar, wind, bind, H_orig, padinfo] = down_recurse(H, b, wbar, Omegas, Qs, wind, bind, H_orig, padinfo)

if H.A11.isleaf
    zk1 = size(H.A12.Z,2);
    zk2 = size(H.A21.Z,2);
    yk1 = size(H.A12.Y,1);
    yk2 = size(H.A21.Y,1);
    m1 = H.A11.size(1);
    m2 = H.A22.size(1);
    n1 = H.A11.size(2);
    n2 = H.A22.size(2);

    % Precondition check: each leaf's own "slack" (its extra width beyond
    % its own row count, n_tau - m_tau) must be large enough to absorb
    % the off-diagonal rank contributed by its SIBLING (A21's rank for
    % A11's slack, A12's rank for A22's slack) -- that's exactly what the
    % Size Reduction LQ step below needs room for. This can fail even for
    % a matrix whose HSS approximation is essentially exact: it is a
    % sizing/geometry precondition (blocksize and aspect ratio vs.
    % off-diagonal rank), not a compression-quality issue.
    if yk2 + m1 > n1
        error('hss_ulvminnormsolve:insufficientSlack', ...
            ['at the leaf pair with rows %s, A11''s off-diagonal rank (from ' ...
            'A21, %d) plus its own row count (%d) exceeds its column count ' ...
            '(%d) -- i.e. A11''s slack (n_tau - m_tau = %d) is smaller than ' ...
            'the rank it needs to absorb. Increase blocksize (or use a ' ...
            'matrix with a wider aspect ratio, or a smaller off-diagonal ' ...
            'rank) so every leaf''s slack exceeds its sibling''s ' ...
            'off-diagonal rank.'], ...
            mat2str(H.A11.Ir), yk2, m1, n1, n1-m1);
    end
    if yk1 + m2 > n2
        error('hss_ulvminnormsolve:insufficientSlack', ...
            ['at the leaf pair with rows %s, A22''s off-diagonal rank (from ' ...
            'A12, %d) plus its own row count (%d) exceeds its column count ' ...
            '(%d) -- i.e. A22''s slack (n_tau - m_tau = %d) is smaller than ' ...
            'the rank it needs to absorb. Increase blocksize (or use a ' ...
            'matrix with a wider aspect ratio, or a smaller off-diagonal ' ...
            'rank) so every leaf''s slack exceeds its sibling''s ' ...
            'off-diagonal rank.'], ...
            mat2str(H.A22.Ir), yk1, m2, n2, n2-m2);
    end

    % Size Reduction
    % [V_tau;D_tau] stacks: tau=1 uses V1=A21.Y (rank yk2) with D1=A11.D
    % (m1 rows); tau=2 uses V2=A12.Y (rank yk1) with D2=A22.D (m2 rows).
    % The kept width after LQ is (own stack's rank + own D's row count).
    [YD1,Omega1s] = lq([H.A12.Y;H.A22.D]);
    [YD2,Omega2s] = lq([H.A21.Y;H.A11.D]);
    H.A12.Y = YD1(1:yk1,1:(yk1+m2));
    H.A21.Y = YD2(1:yk2,1:(yk2+m1));
    H.A22.D = YD1((yk1+1):end,1:(yk1+m2));
    H.A11.D = YD2((yk2+1):end,1:(yk2+m1));
    Omegas = [Omegas,Omega2s,Omega1s];

    % Decoupling

    % QL on Z
    [P1,hatZ1] = ql(H.A12.Z);
    [P2,hatZ2] = ql(H.A21.Z);
    % update Z and D
    H.A12.Z = hatZ1;
    H.A21.Z = hatZ2;
    H.A11.D = P1'*H.A11.D;
    H.A22.D = P2'*H.A22.D;
    % update RHS
    b(1:size(P1,1)) = P1'*b(1:size(P1,1));
    % b(m1+1:size(P2,1)) = P2'*b(m1+1:size(P2,1));
    b(m1+1:end) = P2'*b(m1+1:end);

    % LQ on top of updated D
    upD1 = H.A11.D(1:(m1-zk1),:);
    upD2 = H.A22.D(1:(m2-zk2),:);
    [D1,Q1s] = lq(upD1);
    [D2,Q2s] = lq(upD2);

    % update D
    H.A11.D(1:(m1-zk1),:) = D1;
    H.A22.D(1:(m2-zk2),:) = D2;
    H.A11.D((m1-zk1+1):end,:) = H.A11.D((m1-zk1+1):end,:) * Q1s';
    H.A22.D((m2-zk2+1):end,:) = H.A22.D((m2-zk2+1):end,:) * Q2s';
    % Q1s comes from tau=1's own D-tilde, so it applies to tau=1's own V
    % (V1 = A21.Y); Q2s applies to tau=2's own V (V2 = A12.Y).
    H.A21.Y = H.A21.Y*Q1s';
    H.A12.Y = H.A12.Y*Q2s';

    Qs = [Qs,Q1s,Q2s];

    % refresh size/Ic bookkeeping now that the column widths have shrunk
    % (needed so H*wbar via hss_matvec stays consistent, and so the
    % ancestor levels read the correct reduced widths for their own
    % index shifting).
    H.A11.size = size(H.A11.D);
    H.A22.size = size(H.A22.D);
    H.A12.size = [size(H.A12.Z,1), size(H.A12.Y,2)];
    H.A21.size = [size(H.A21.Z,1), size(H.A21.Y,2)];
    H.A11.Ic = [H.Ic(1), H.Ic(1)+H.A11.size(2)-1];
    H.A22.Ic = [H.A11.Ic(2)+1, H.A11.Ic(2)+H.A22.size(2)];
    H.A12.Ic = H.A22.Ic;  H.A12.Ir = H.A11.Ir;
    H.A21.Ic = H.A11.Ic;  H.A21.Ir = H.A22.Ir;
    H.Ic = [H.A11.Ic(1), H.A22.Ic(2)];
    H.size = [H.A11.size(1)+H.A22.size(1), H.A11.size(2)+H.A22.size(2)];

    redw1 = H.A11.size(2);
    redw2 = H.A22.size(2);

    % Solving -- wbar from here on lives in the size-reduced coordinate
    % system (total length redw1+redw2), not the original width, so its
    % offsets use redw1 rather than the pre-reduction n1.
    wbar = zeros(redw1+redw2,1);
    wbar(1:(m1-zk1)) = H.A11.D(1:(m1-zk1),1:(m1-zk1))\b(1:(m1-zk1));
    wbar((redw1+1):(redw1+m2-zk2)) = H.A22.D(1:(m2-zk2),1:(m2-zk2))\b((m1+1):(m1+m2-zk2));
    bind1 = 1:(m1-zk1);
    bind2 = (m1+1):(m1+m2-zk2);
    bind = [bind1,bind2];
    wind1 = 1:(m1-zk1);
    wind2 = (redw1+1):(redw1+m2-zk2);
    wind = [wind1,wind2];

    % track (original width, reduced width) per leaf, in the same order
    % Omegas/Qs are collected, so the outer function can pad-then-lift
    % each leaf's chunk back through its own Omega at the end.
    padinfo = [padinfo; n1, redw1; n2, redw2];

else
    m1 = H.A11.size(1);
    n1 = H.A11.size(2);
    [H.A11,Omegas,Qs,brow1,wbarseg1,wind1,bind1,H_orig.A11,padinfo] = down_recurse(H.A11, b(1:m1), wbar(1:n1), Omegas, Qs, [],[],H_orig.A11,padinfo);
    [H.A22,Omegas,Qs,brow2,wbarseg2,wind2,bind2,H_orig.A22,padinfo] = down_recurse(H.A22, b(m1+1:end), wbar(n1+1:end), Omegas, Qs, [],[],H_orig.A22,padinfo);
    wshift = H.A11.size(2);
    bshift = H.A11.size(1);
    b = [brow1; brow2];
    wbar = [wbarseg1; wbarseg2];
    wind = [wind1,wind2+wshift];
    bind = [bind1,bind2+bshift];
    H.size = [H.A11.size(1)+H.A22.size(1), H.A11.size(2)+H.A22.size(2)];
    H.Ic = [H.A11.Ic(1), H.A22.Ic(2)];
    % off-diagonal Ic/size at this level: A12 spans rows of A11, cols of A22
    H.A12.Ic = H.A22.Ic;  H.A12.Ir = H.A11.Ir;
    H.A21.Ic = H.A11.Ic;  H.A21.Ir = H.A22.Ir;
    H.A12.size = [H.A11.size(1), H.A22.size(2)];
    H.A21.size = [H.A22.size(1), H.A11.size(2)];
    % H.Ir is unchanged
end
end


function tf = slack_ok(H)
% Same precondition down_recurse's leaf-pair branch guards against
% (yk2+m1<=n1 and yk1+m2<=n2), checked recursively over the whole tree so
% the caller can decide, up front, whether the ULV reduction can run to
% completion without hitting the insufficient-slack case anywhere.
if H.isleaf
    tf = true;
    return
end
if H.A11.isleaf
    yk1 = size(H.A12.Y,1);
    yk2 = size(H.A21.Y,1);
    m1 = H.A11.size(1);
    m2 = H.A22.size(1);
    n1 = H.A11.size(2);
    n2 = H.A22.size(2);
    tf = (yk2 + m1 <= n1) && (yk1 + m2 <= n2);
else
    tf = slack_ok(H.A11) && slack_ok(H.A22);
end
end


function [Q,H] = leafcompression(H)
% identical to hss_ulvvecsolve.m's leafcompression: QL of the off-diagonal
% row generator Z, via flip+QR. Z = QRF = (QF)(FRF) -> (QF)'Z = FRF = newZ
Fleft = flip(eye(size(H.Z,1)));
Fright = flip(eye(size(H.Z,2)));
[Q,R] = qr(H.Z*Fright);
H.Z = Fleft * (R) * Fright;
Q = (Q*Fleft);
end


function [Qd,H] = rowdecouple(H,k)
% LQ of the top (n-k) rows of D (own column block only), applied as a
% right-multiply to the FULL D. The top rows of the result are then
% [triangular, 0] -- decoupled from the free (last 2k) columns -- because
% L from an LQ of a wide block is lower triangular and so has zero
% columns past its row count.
n = H.size(1);
top = H.D(1:n-k,:);
[~,Qd] = lq(top);
H.D = H.D*Qd';
end


function [wbar,wind,rowind] = bsolve_minnorm(D,b,k)
% D is n x (n+k) here (n rows, k off-diagonal rank, so n+k columns after
% size reduction). Top (n-k) rows/cols are triangular and uniquely solve;
% bottom k rows and last 2k columns remain free for the level above.
[n,kn] = size(D);
t = n-k;
w = D(1:t,1:t)\b(1:t);
wbar = zeros(kn,1);
wind = false(kn,1);
wbar(1:t) = w;
wind(1:t) = true;
rowind = false(n,1);
rowind(t+1:end) = true;
end


function H = discard(H)
if H.A11.isleaf
    [m1,n1] = size(H.A11.D);
    [m2,n2] = size(H.A22.D);
    k1 = size(H.A12.Z,2);
    k2 = size(H.A21.Z,2);

    H.A12.Z = H.A12.Z(end-k1+1:end,:);
    H.A12.Y = H.A12.Y(:,m2-k2+1:end);
    H.A12.size = [size(H.A12.Z,1), size(H.A12.Y,2)];
    H.A12.isleaf = 1;

    H.A21.Z = H.A21.Z(end-k2+1:end,:);
    H.A21.Y = H.A21.Y(:,m1-k1+1:end);
    H.A21.size = [size(H.A21.Z,1), size(H.A21.Y,2)];
    H.A21.isleaf = 1;

    H.A11.D = H.A11.D(m1-k1+1:end, m1-k1+1:end);
    H.A11.size = size(H.A11.D);
    H.A22.D = H.A22.D(m2-k2+1:end, m2-k2+1:end);
    H.A22.size = size(H.A22.D);
else
    H.A11 = discard(H.A11);
    H.A22 = discard(H.A22);
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
    H.size = size(H.D);
    H.levelcount = H.levelcount-1;
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
% UNCHANGED from hss_ulvvecsolve.m.
if H.isdiag
    H.D = [H.A11.D, H.A12.Z*H.A12.lrcomponent*H.A12.Y; H.A21.Z*H.A21.lrcomponent*H.A21.Y,H.A22.D];
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
    Z1 = Hparent.A11.A12.Z;
    Z2 = Hparent.A11.A21.Z;
    Hparent.A12.Z = blkdiag(Z1,Z2)*Hparent.A12.Z;

    Y1 = Hparent.A22.A21.Y;
    Y2 = Hparent.A22.A12.Y;
    Hparent.A12.Y = Hparent.A12.Y*blkdiag(Y1,Y2);
    H.levelcount = H.levelcount-1;
    Hparent.A12.size = [size(Hparent.A12.Z,1),size(Hparent.A12.Y,2)];

elseif mod(H.rowtreeindex,2) == 0
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


function [Q, L] = ql(A)
[Q_r, R] = qr(flipud(fliplr(A)));
L = flipud(fliplr(R));
Q = flipud(fliplr(Q_r));
end

function [L, Q] = lq(A)
[Q_t, R] = qr(A');
L = R';
Q = Q_t';
end