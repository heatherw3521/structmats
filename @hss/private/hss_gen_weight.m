function G = hss_gen_weight(G, S, L)
%HSS_GEN_WEIGHT  Generators of S * H * inv(L) for leaf-conforming S and L.
%   S (rows) and L (columns) are each [] (identity), a vector (diagonal), or a
%   cell array of square blocks, one per leaf, in leaf order. Only the leaf
%   generators change (memo 2, Lemma 3.1):
%       D_tau <- S_tau * D_tau / L_tau,   U_tau <- S_tau * U_tau,
%       V_tau <- L_tau' \ V_tau          (V_tau^* <- V_tau^* / L_tau).
%   L is applied by solves, never by forming inv(L_tau).
nl = numel(G.D);
rb = G.rb{G.L+1}; cb = G.cb{G.L+1};
for i = 1:nl
    r = rb(i)+1:rb(i+1);  c = cb(i)+1:cb(i+1);
    if ~isempty(S)
        if iscell(S)
            G.D{i} = S{i} * G.D{i};
            if G.L > 0, G.U{G.L+1}{i} = S{i} * G.U{G.L+1}{i}; end
        else
            G.D{i} = S(r(:)) .* G.D{i};
            if G.L > 0, G.U{G.L+1}{i} = S(r(:)) .* G.U{G.L+1}{i}; end
        end
    end
    if ~isempty(L)
        if iscell(L)
            G.D{i} = G.D{i} / L{i};
            if G.L > 0, G.V{G.L+1}{i} = L{i}' \ G.V{G.L+1}{i}; end
        else
            G.D{i} = G.D{i} ./ reshape(L(c), 1, []);
            if G.L > 0, G.V{G.L+1}{i} = G.V{G.L+1}{i} ./ conj(L(c(:))); end
        end
    end
end
end
