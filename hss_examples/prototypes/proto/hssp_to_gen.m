function G = hssp_to_gen(H)
%HSSP_TO_GEN  Read the generators (D, U/R, V/W, B12, B21, cluster boundaries)
%   out of an hss object; inverse of hss_from_gen. Requires every leaf at the same depth (true for
%   everything hss_constructor.m builds); otherwise errors.
L = H.levelcount;
G.L = L;
G.rb = cell(L+1,1); G.cb = cell(L+1,1);
G.U = cell(L+1,1); G.V = cell(L+1,1);
G.B12 = cell(L,1); G.B21 = cell(L,1);
G.D = cell(2^L,1);
for l = 0:L
    G.rb{l+1} = zeros(1, 2^l+1); G.cb{l+1} = zeros(1, 2^l+1);
    G.U{l+1} = cell(2^l,1); G.V{l+1} = cell(2^l,1);
    if l < L, G.B12{l+1} = cell(2^l,1); G.B21{l+1} = cell(2^l,1); end
end
if L == 0
    G.D{1} = H.D; G.rb{1} = [0 size(H, 1)]; G.cb{1} = [0 size(H, 2)];
    return
end
G = walk(H, 0, 1, G);
if any(cellfun(@isempty, G.D))
    error('hss:nonuniformTree', ['weighted minnorm/tikhonov need every leaf at the ' ...
        'same depth (levelcount = %d); this tree is not uniform.'], L);
end
end

function G = walk(N, l, i, G)
G.rb{l+1}(i) = N.Ir(1)-1;  G.rb{l+1}(i+1) = N.Ir(2);
G.cb{l+1}(i) = N.Ic(1)-1;  G.cb{l+1}(i+1) = N.Ic(2);
if N.isleaf
    if l == G.L, G.D{i} = N.D; end     % a leaf above level L leaves a gap (checked)
    return
end
G.B12{l+1}{i} = N.A12.lrcomponent;
G.B21{l+1}{i} = N.A21.lrcomponent;
G.U{l+2}{2*i-1} = N.A12.Z;   G.V{l+2}{2*i}   = N.A12.Y';
G.U{l+2}{2*i}   = N.A21.Z;   G.V{l+2}{2*i-1} = N.A21.Y';
G = walk(N.A11, l+1, 2*i-1, G);
G = walk(N.A22, l+1, 2*i, G);
end
