function S = mp_exp_sweeps(opts)
%MP_EXP_SWEEPS  Time of H\b vs HSS rank, vs blocksize, and vs aspect ratio
%   m/n (Assumption 1 of memo 1 fails at large m/n; the solver does not need
%   it).  F2 exact random HSS.  Not used in memo 1.
%   Times are of the first H\b (factor + solve), see mp_time_solve.
if nargin < 1, opts = struct(); end
q = isfield(opts,'quick') && opts.quick;
reps = 3; if q, reps = 1; end
rng(5);
n = 2^16; if q, n = 2^12; end
S.rank = [];
% two FIXED leaf sizes: blocksize 64 for k <= 32, blocksize 128 for k <= 64
for cfg = {struct('bs', 64, 'ks', [2 4 8 16 24 32]), struct('bs', 128, 'ks', [2 4 8 16 32 48 64])}
    bs = cfg{1}.bs;
    for k = cfg{1}.ks
        G = mp_gen_random(n/2, n, k, bs); H = mp_hss_from_generators(G);
        [b, xs] = mp_manufactured(H);
        [t, ~, x] = mp_time_solve(H, b, reps);       % first solve (factor + solve)
        S.rank = [S.rank, struct('k', k, 'bs', bs, 'L', G.L, 't', t, 'err', mp_rel(x, xs))];
        fprintf('k=%3d bs=%3d  H\\b %.3fs  err %.1e\n', k, bs, t, mp_rel(x, xs));
    end
end
S.blocksize = [];
for bs = [16 32 64 128 256 512]
    G = mp_gen_random(n/2, n, 10, bs); H = mp_hss_from_generators(G);
    [b, xs] = mp_manufactured(H);
    [t, ~, x] = mp_time_solve(H, b, reps);       % first solve (factor + solve)
    S.blocksize = [S.blocksize, struct('bs', bs, 'L', G.L, 't', t, 'err', mp_rel(x, xs))];
    fprintf('bs=%3d L=%2d  H\\b %.3fs  err %.1e\n', bs, G.L, t, mp_rel(x, xs));
end
n = 2^14; if q, n = 2^11; end
S.aspect = [];
for ratio = [0.1 0.25 0.5 0.6 0.7 0.75 0.8 0.85 0.9 0.95]
    m = round(ratio*n);
    G = mp_gen_random(m, n, 16, 64); H = mp_hss_from_generators(G);
    [b, xs] = mp_manufactured(H);
    [t, ~, x] = mp_time_solve(H, b, reps);       % first solve (factor + solve)
    [ok, nbad] = mp_slack_ok(G);
    S.aspect = [S.aspect, struct('ratio', ratio, 'm', m, 'n', n, 't', t, 'err', mp_rel(x, xs), ...
                                 'slack_ok', ok, 'leaves_violating', nbad)];
    fprintf('m/n=%.2f  slack_ok=%d (%d leaves violate)  H\\b %.3fs  err %.1e\n', ratio, ok, nbad, t, mp_rel(x, xs));
end
end
