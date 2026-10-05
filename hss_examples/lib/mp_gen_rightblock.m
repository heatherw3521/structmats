function G = mp_gen_rightblock(G, M)
%MP_GEN_RIGHTBLOCK  Generators of H * blkdiag(M{1},...,M{2^L}) for square
%   blocks M{i} conforming to the leaf column partition:
%   D_tau <- D_tau M_tau,  V_tau <- M_tau' V_tau.
L = G.L;
for i = 1:2^L
    G.D{i} = G.D{i} * M{i};
    if L > 0, G.V{L+1}{i} = M{i}' * G.V{L+1}{i}; end
end
end
