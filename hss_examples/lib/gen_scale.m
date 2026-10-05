function G = gen_scale(G, s_rows, s_cols)
%GEN_SCALE  Generators of diag(s_rows) * H * diag(s_cols) (either may be []).
%   Only the leaf generators change:  D_tau <- S_tau D_tau C_tau,
%   U_tau <- S_tau U_tau,  V_tau <- conj(C_tau) V_tau.  Same tree, same ranks.
L = G.L;
for i = 1:2^L
    r = G.rb{L+1}(i)+1:G.rb{L+1}(i+1); c = G.cb{L+1}(i)+1:G.cb{L+1}(i+1);
    if ~isempty(s_rows)
        G.D{i} = s_rows(r(:)) .* G.D{i};
        if L > 0, G.U{L+1}{i} = s_rows(r(:)) .* G.U{L+1}{i}; end
    end
    if ~isempty(s_cols)
        G.D{i} = G.D{i} .* reshape(s_cols(c), 1, []);
        if L > 0, G.V{L+1}{i} = conj(s_cols(c(:))) .* G.V{L+1}{i}; end
    end
end
end
