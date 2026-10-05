function [ok, nbad] = mp_slack_ok(G)
%MP_SLACK_OK  Slack condition on the leaf level of generators G:
%   l_tau + n_tau <= p_tau for every leaf (l_tau = column rank, n_tau rows,
%   p_tau columns). Reported by the experiments for reference only; the
%   solver does not need it.
L = G.L; nbad = 0;
if L == 0, ok = true; return, end
for i = 1:2^L
    [n, p] = size(G.D{i});
    if size(G.V{L+1}{i}, 2) + n > p, nbad = nbad + 1; end
end
ok = (nbad == 0);
end
