function H = hss_from_gen(G)
%HSS_FROM_GEN  Assemble an @hss object directly from HSS generators,
%   never forming the dense matrix.  The object has exactly the layout
%   hss_constructor.m produces (same fields, tree indices and nested Z/Y/
%   lrcomponent convention), so H*x, H'*x, H\b, full(H) all work unchanged.
%
%   G.L            leaf level (>= 1)
%   G.rb{l+1}      0-based row boundaries of the level-l clusters (length 2^l+1)
%   G.cb{l+1}      0-based column boundaries
%   G.D{i}         leaf diagonal block i (i = 1..2^L)
%   G.U{l+1}{i}    l = L : leaf row basis (n_i x k_i);  1 <= l < L : translation R
%   G.V{l+1}{i}    l = L : leaf column basis (p_i x l_i); 1 <= l < L : translation W
%   G.B12{l+1}{i}, G.B21{l+1}{i}  (l = 0..L-1) couplings between the children of
%                  node (l,i):  H(I_a,J_c) = Ubig_a*B12*Vbig_c',  H(I_c,J_a) = Ubig_c*B21*Vbig_a'
%   Optional: G.blocksize (stored at the root, informational).
%
%   Correspondence with the class:  Z = U (or R),  Y = V' (or W'),  lrcomponent = B.
%   (Same code as hss_examples/minnorm_paper/lib/mp_hss_from_generators.m.)
L = G.L;
m = G.rb{1}(end); n = G.cb{1}(end);
H = hss();
H.sz = [m n];
H.isroot = true;
H.isdiag = true;
H.level = 0;
H.rowtreeindex = 1;
H.coltreeindex = 1;
H.Ir = [1 m];
H.Ic = [1 n];
H.levelcount = L;
if isfield(G,'blocksize'), H.blocksize = G.blocksize; end
if L == 0
    H.isleaf = 1;
    H.D = G.D{1};
    return
end
H.isleaf = false;
H.A11 = diagnode(G, 1, 1);
H.A22 = diagnode(G, 1, 2);
H.A12 = offnode(G, 1, 1, 2, G.B12{1}{1});
H.A21 = offnode(G, 1, 2, 1, G.B21{1}{1});
end

function N = diagnode(G, l, i)
L = G.L;
N = hss();
N.isroot = false;
N.isdiag = true;
N.level = l;
N.rowtreeindex = i;
N.coltreeindex = i;
N.Ir = [G.rb{l+1}(i)+1, G.rb{l+1}(i+1)];
N.Ic = [G.cb{l+1}(i)+1, G.cb{l+1}(i+1)];
N.sz = [N.Ir(2)-N.Ir(1)+1, N.Ic(2)-N.Ic(1)+1];
N.levelcount = L;
if l == L
    N.isleaf = true;
    N.D = G.D{i};
else
    N.isleaf = false;
    N.A11 = diagnode(G, l+1, 2*i-1);
    N.A22 = diagnode(G, l+1, 2*i);
    N.A12 = offnode(G, l+1, 2*i-1, 2*i, G.B12{l+1}{i});
    N.A21 = offnode(G, l+1, 2*i, 2*i-1, G.B21{l+1}{i});
end
end

function N = offnode(G, l, ir, ic, B)
% off-diagonal block at level l: rows of cluster (l,ir), columns of (l,ic)
N = hss();
N.isroot = false;
N.isdiag = false;
N.level = l;
N.rowtreeindex = ir;
N.coltreeindex = ic;
N.Ir = [G.rb{l+1}(ir)+1, G.rb{l+1}(ir+1)];
N.Ic = [G.cb{l+1}(ic)+1, G.cb{l+1}(ic+1)];
N.sz = [N.Ir(2)-N.Ir(1)+1, N.Ic(2)-N.Ic(1)+1];
N.isleaf = (l == G.L);
N.Z = G.U{l+1}{ir};
N.Y = G.V{l+1}{ic}';
N.lrcomponent = B;
N.lowrankrows = [];
N.lowrankcols = [];
end
