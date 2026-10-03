function [Ga, ix, is] = hss_gen_augment(G, lam)
%HSS_GEN_AUGMENT  Generators of the Tikhonov-augmented matrix [H, lam*I], with
%   the identity columns belonging to leaf tau's rows placed right after leaf
%   tau's own columns (a column permutation).  The result is HSS on the SAME
%   tree with the SAME ranks:  D_tau -> [D_tau, lam*I],  V_tau -> [V_tau; 0].
%   ix / is : positions of the x-part and of the s-part in the augmented vector.
L = G.L;
Ga = G;
nl = 2^L;
cbL = zeros(1, nl+1); ix = []; is = [];
off = 0;
for i = 1:nl
    ni = G.rb{L+1}(i+1)-G.rb{L+1}(i); pi_ = G.cb{L+1}(i+1)-G.cb{L+1}(i);
    Ga.D{i} = [G.D{i}, lam*eye(ni)];
    if L > 0
        Ga.V{L+1}{i} = [G.V{L+1}{i}; zeros(ni, size(G.V{L+1}{i},2))];
    end
    ix = [ix, off+(1:pi_)]; %#ok<AGROW>
    is = [is, off+pi_+(1:ni)]; %#ok<AGROW>
    off = off + pi_ + ni;
    cbL(i+1) = off;
end
Ga.cb{L+1} = cbL;
for l = L-1:-1:0
    Ga.cb{l+1} = Ga.cb{l+2}(1:2:end);
end
end
