function H = hss_matmat(A,B)
%HSS_MATMAT product of two HSS matrices, A*B, returned as a new HSS matrix.
%
% Works at any depth. The hard part beyond levelcount<=1 is that
% off-diagonal generators above the leaf-pair level are nested/telescoping
% TRANSLATION matrices (mapping a child's compressed basis up to its
% parent's), not direct compressions of the dense off-diagonal block --
% see hss_matvec's ascend/descend for the same structure applied to a
% vector. Multiplying two such matrices while staying in this nested form
% (never materializing an O(N)-sized basis) needs:
%
%   1. Gamma(p) := V_p^full * U_p^full, a small cross-Gram matrix between
%      A's column-basis and B's row-basis at every diagonal node p,
%      computed bottom-up (gammapair below) purely from each level's own
%      translation/lrcomponent fields -- no dependence on the product
%      tree being built at all.
%   2. The diagonal recursion C.p = A.p*B.p + correction, where
%      "correction" is a low-rank term living outside p's own index range
%      (A12*B21-type cross terms from every ancestor level). Rather than
%      re-deriving that correction from scratch at each level, it is
%      threaded down as a small core matrix Gin together with the two
%      external translations (ZextA, YextB) that embed p's own basis into
%      the frame where Gin applies -- see matmatrec.
%   3. The off-diagonal construction (buildoffdiag), which generalizes the
%      levelcount==1 formula below: C's rank at every off-diagonal is
%      simply A's rank plus B's rank (concatenation, no truncation), and
%      when the off-diagonal is itself non-leaf, both the "direct" half
%      (A's or B's own translation) and the "crossed" half (derived from
%      Gamma at the children) have to be re-expressed in the interleaved
%      [A-part,B-part] basis that the recursion below builds for every
%      node -- embedrows/embedcols do that zero-padding.
%
% Setting Gin=[] throughout (the top-level call) makes step 2 a no-op and
% the whole recursion collapses exactly to the original levelcount==1
% formula at every leaf-pair parent -- that formula is kept below as the
% base case docstring for reference. Verified against a dense reference
% (real and complex, square and rectangular, levelcount up to several
% levels deep) via a from-scratch derivation checked against random
% synthetic HSS trees before being ported here; see hss_examples/tests for
% the corresponding regression tests.

if A.size(2) ~= B.size(1)
    error("hss_matmat:sizeMismatch", "Mismatched sizes for matrix multiplication.")
elseif A.levelcount ~= B.levelcount
    error("hss_matmat:levelcountMismatch", ...
        "HSS matrices must have the same structure for efficient matmat multiplication.")
end
checkconformal(A,B);

H = matmatrec(A,B,[],[],[]);

end

% =========================================================================
function checkconformal(A,B)
% Matching top-level size and levelcount are NOT enough on their own: two
% HSS trees built with different blocksize/cutrule choices can share both
% while splitting their shared dimension (A's columns / B's rows) at
% different points, which would otherwise silently misalign every
% recursive call below (A.A11 multiplying against a B.A11 that doesn't
% actually correspond to the same index range) rather than erroring
% cleanly. Walk both trees together and check every split lines up.

if A.isleaf ~= B.isleaf
    error('hss_matmat:incompatibleStructure', ...
        ['A and B have incompatible HSS tree structures: one is a leaf ' ...
        'at a level where the other is not, even though both report the ' ...
        'same levelcount.']);
end
if A.isleaf
    return
end
if A.A11.size(2) ~= B.A11.size(1)
    error('hss_matmat:incompatibleStructure', ...
        ['A and B have incompatible HSS tree structures: A''s column ' ...
        'partition does not match B''s row partition along the shared ' ...
        'dimension (A.A11 has %d columns, B.A11 has %d rows). Both ' ...
        'matrices must be built with the same split points -- e.g. the ' ...
        'same blocksize and cutrule -- for their product to have a ' ...
        'well-defined HSS structure.'], A.A11.size(2), B.A11.size(1));
end
checkconformal(A.A11,B.A11);
checkconformal(A.A22,B.A22);

end

% =========================================================================
function H = matmatrec(A,B,Gin,ZextA,YextB)
% Core recursive HSS x HSS product. Gin (or [] if none) is a small core
% matrix such that A/B's parent wants an extra correction
% ZextA*Gin*YextB added to A*B, where ZextA/YextB are the translations
% embedding A's/B's own basis into the ancestor frame where Gin lives (see
% header comment above). At a leaf this correction is applied directly;
% at a non-leaf it is split across both diagonal children (folded together
% with each node's own native cross term, from Gamma) and contributes a
% single small cross-block to the off-diagonal lrcomponent.

if A.isleaf
    H = A;
    H.Ic = B.Ic;   % rows follow A's own structure (kept via H=A); columns follow B's
    H.D = A.D*B.D;
    if ~isempty(Gin)
        H.D = H.D + ZextA*Gin*YextB;
    end
    H.size = size(H.D);
    return
end

c1A = A.A11; c2A = A.A22;
c1B = B.A11; c2B = B.A22;

% Gamma at this level's two children (bottom-up; g1/g2 are [] for a leaf
% child, else the {Gamma(child.A11),Gamma(child.A22)} pair one level
% further down, needed again below to build the off-diagonal nodes).
if c1A.isleaf
    g1 = [];
    Gc1 = A.A21.Y * B.A12.Z;
else
    [g1a,g1b] = gammapair(c1A,c1B);
    g1 = {g1a,g1b};
    Gc1 = A.A21.Y * blkdiag(g1a,g1b) * B.A12.Z;
end
if c2A.isleaf
    g2 = [];
    Gc2 = A.A12.Y * B.A21.Z;
else
    [g2a,g2b] = gammapair(c2A,c2B);
    g2 = {g2a,g2b};
    Gc2 = A.A12.Y * blkdiag(g2a,g2b) * B.A21.Z;
end

% native cross-term correction (A12*B21-type) generated at THIS level
Gnative_c1 = A.A12.lrcomponent * Gc2 * B.A21.lrcomponent;
Gnative_c2 = A.A21.lrcomponent * Gc1 * B.A12.lrcomponent;

cross_c1c2 = [];
cross_c2c1 = [];
if isempty(Gin)
    Gtotal_c1 = Gnative_c1;
    Gtotal_c2 = Gnative_c2;
else
    kAc1 = size(A.A12.Z,2);
    kBc1 = size(B.A21.Y,1);
    ZextA_top = ZextA(1:kAc1,:);       ZextA_bot = ZextA(kAc1+1:end,:);
    YextB_left = YextB(:,1:kBc1);      YextB_right = YextB(:,kBc1+1:end);
    Gtotal_c1 = Gnative_c1 + ZextA_top*Gin*YextB_left;
    Gtotal_c2 = Gnative_c2 + ZextA_bot*Gin*YextB_right;
    cross_c1c2 = ZextA_top*Gin*YextB_right;
    cross_c2c1 = ZextA_bot*Gin*YextB_left;
end

H_c1 = matmatrec(c1A,c1B,Gtotal_c1,A.A12.Z,B.A21.Y);
H_c2 = matmatrec(c2A,c2B,Gtotal_c2,A.A21.Z,B.A12.Y);

H12 = buildoffdiag(c1A,c1B,c2A,c2B,A.A12,B.A12,g1,g2,cross_c1c2);
H21 = buildoffdiag(c2A,c2B,c1A,c1B,A.A21,B.A21,g2,g1,cross_c2c1);

% H_c1/H_c2/H12/H21 already carry correct GLOBAL Ir (inherited from A's
% own tree, unchanged) and Ic (inherited from B's own tree -- see the leaf
% branch above and buildoffdiag below) -- no re-derivation needed, and
% re-deriving them "from 1" here would be wrong for anything but the very
% top-level call (this node's own Ir/Ic are only [1,...] when A/B are the
% true root).
H = A;
H.Ic = B.Ic;
H.A11 = H_c1;
H.A22 = H_c2;
H.A12 = H12;
H.A21 = H21;

H.size = [H.A11.size(1)+H.A22.size(1), H.A11.size(2)+H.A22.size(2)];

end

% =========================================================================
function [g1,g2] = gammapair(A,B)
% Gamma(A.A11,B.A11), Gamma(A.A22,B.A22) -- the cross-Gram
% V_p^full(A-tree) * U_p^full(B-tree) at A's two children, computed
% bottom-up. A,B must be matching non-leaf diagonal nodes (same tree
% shape). Depends only on A's and B's own stored fields, never on the
% product tree under construction.

if A.A11.isleaf
    g1 = A.A21.Y * B.A12.Z;
else
    [g1a,g1b] = gammapair(A.A11,B.A11);
    g1 = A.A21.Y * blkdiag(g1a,g1b) * B.A12.Z;
end
if A.A22.isleaf
    g2 = A.A12.Y * B.A21.Z;
else
    [g2a,g2b] = gammapair(A.A22,B.A22);
    g2 = A.A12.Y * blkdiag(g2a,g2b) * B.A21.Z;
end

end

% =========================================================================
function off = buildoffdiag(c1A,c1B,c2A,c2B,Aoff,Boff,g1pair,g2pair,crossterm)
% The product's off-diagonal node pairing rows=c1,cols=c2 (i.e. C.A12 if
% c1=A.A11, or C.A21 if c1=A.A22 -- called once for each). Aoff/Boff are
% the corresponding A-tree/B-tree offdiag nodes (parent.A12, say); g1pair
% = {Gamma(c1.A11),Gamma(c1.A22)} (or [] if c1 is a leaf), g2pair likewise
% for c2 (both already computed by the caller). crossterm is the (c1,c2)
% block of an inherited correction from an ancestor level, or [] if none.
%
% Rank adds (no truncation): Z = [direct-A part, B-crossed part], Y =
% [A-crossed part; direct-B part], lrcomponent = [[A's lr, crossterm];
% [0, B's lr]] -- this is exactly the levelcount==1 formula in hss_matmat's
% header when c1/c2 are leaves (direct = A's/B's own leaf Z/Y, crossed =
% A11.D*B12.Z / A12.Y*B22.D, no interleaving needed since a leaf has no
% children to interleave). When c1 or c2 is itself non-leaf, both halves
% instead have to be produced in the SAME interleaved [A-part,B-part]
% child ordering used everywhere else in this file (embedrows/embedcols
% zero-pad the "direct" half into it; the "crossed" half comes out already
% in that order from the Gamma-based recursion).

if c1A.isleaf
    Z1 = c1A.D * Boff.Z;
    Zdirect = Aoff.Z;
else
    g1a = g1pair{1}; g1b = g1pair{2};
    kAc1c1 = size(c1A.A12.Z,2); kAc1c2 = size(c1A.A21.Z,2);
    kBc1c1 = size(c1B.A12.Z,2); kBc1c2 = size(c1B.A21.Z,2);
    RBc1 = Boff.Z;
    RBc1_top = RBc1(1:kBc1c1,:);
    RBc1_bot = RBc1(kBc1c1+1:end,:);
    block_c1c1 = [c1A.A12.lrcomponent*g1b*RBc1_bot; RBc1_top];
    block_c1c2 = [c1A.A21.lrcomponent*g1a*RBc1_top; RBc1_bot];
    Z1 = [block_c1c1; block_c1c2];
    Zdirect = embedrows(Aoff.Z, [kAc1c1,kAc1c2], [kBc1c1,kBc1c2]);
end

if c2A.isleaf
    Y2 = Aoff.Y * c2B.D;
    Ydirect = Boff.Y;
else
    g2a = g2pair{1}; g2b = g2pair{2};
    kAc2c1 = size(c2A.A21.Y,1); kAc2c2 = size(c2A.A12.Y,1);
    kBc2c1 = size(c2B.A21.Y,1); kBc2c2 = size(c2B.A12.Y,1);
    RAc2 = Aoff.Y;
    RAc2_left = RAc2(:,1:kAc2c1);
    RAc2_right = RAc2(:,kAc2c1+1:end);
    coef1 = [RAc2_left, RAc2_right*g2b*c2B.A21.lrcomponent];
    coef2 = [RAc2_right, RAc2_left*g2a*c2B.A12.lrcomponent];
    Y2 = [coef1, coef2];
    Ydirect = embedcols(Boff.Y, [kAc2c1,kAc2c2], [kBc2c1,kBc2c2]);
end

off = Aoff;
off.Ic = Boff.Ic;   % rows (Ir) follow A's own structure; columns follow B's
off.Z = [Zdirect, Z1];
off.Y = [Y2; Ydirect];

kA1 = size(Aoff.lrcomponent,1); kA2 = size(Aoff.lrcomponent,2);
kB1 = size(Boff.lrcomponent,1); kB2 = size(Boff.lrcomponent,2);
if isempty(crossterm)
    crossterm = zeros(kA1,kB2);
end
off.lrcomponent = [Aoff.lrcomponent, crossterm; zeros(kB1,kA2), Boff.lrcomponent];
off.isleaf = c1A.isleaf;
off.size = [c1A.size(1), c2B.size(2)];

end

% =========================================================================
function M = embedrows(M0, kAlist, kBlist)
% Zero-pad/interleave M0 (rows grouped only by kAlist, one group per
% child) into the combined [A-part;B-part]-per-child row space sized by
% kAlist+kBlist -- i.e. place each of M0's row-groups into the A-slot of
% its child, zero elsewhere.
total = sum(kAlist) + sum(kBlist);
if isreal(M0)
    M = zeros(total, size(M0,2));
else
    M = complex(zeros(total, size(M0,2)));
end
rout = 0; rin = 0;
for i = 1:numel(kAlist)
    a = kAlist(i); b = kBlist(i);
    M(rout+1:rout+a,:) = M0(rin+1:rin+a,:);
    rout = rout+a+b;
    rin = rin+a;
end
end

% =========================================================================
function M = embedcols(M0, kAlist, kBlist)
% As embedrows, but for columns, and M0's groups land in the B-slot
% (second within each child's [A-part,B-part] column range).
total = sum(kAlist) + sum(kBlist);
if isreal(M0)
    M = zeros(size(M0,1), total);
else
    M = complex(zeros(size(M0,1), total));
end
cout = 0; cin = 0;
for i = 1:numel(kAlist)
    a = kAlist(i); b = kBlist(i);
    cout = cout+a;
    M(:,cout+1:cout+b) = M0(:,cin+1:cin+b);
    cout = cout+b;
    cin = cin+b;
end
end
