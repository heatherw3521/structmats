function S = mp_exp_scaling(opts)
%MP_EXP_SCALING  E3/E4: time and accuracy of
%   H\b as n grows, on manufactured problems whose minimum-norm solution is
%   known (x* = H'z, b = H x*), so no dense matrix is ever formed at large n.
%   For F2 it also times the dense QR minimum-norm solve and CGNE (pcg on
%   H*H' with HSS products) on the same problems.
%
%   Each H\b is timed twice (mp_time_solve): the first solve, which factors
%   and solves, and a repeat solve, which uses the factors H keeps.
%   Leaves use the hss constructor's default blocksize (200) and its depth
%   rule (mp_tree).
%
%   opts.maxexp   (default 17): largest n = 2^maxexp columns (memory grows
%                 like n * blocksize; 2^17 needs a few GB)
%   opts.maxdense (default 13): largest n for the dense QR comparison
%   opts.maxcg    (default 15): largest n for the CGNE comparison
%   opts.reps     (default 3) : timing repetitions (median), after a warm-up
%   opts.quick    : n up to 2^12, one repetition
%   Writes results/scaling.txt (one line per family, method and n).
if nargin < 1, opts = struct(); end
if ~isfield(opts, 'maxexp'), opts.maxexp = 17; end
if ~isfield(opts, 'maxdense'), opts.maxdense = 13; end
if ~isfield(opts, 'maxcg'), opts.maxcg = 15; end
if ~isfield(opts, 'reps'), opts.reps = 3; end
q = isfield(opts, 'quick') && opts.quick;
if q, opts.maxexp = 12; opts.maxdense = 11; opts.maxcg = 11; opts.reps = 1; end
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
bs = 200;                                   % hss_constructor default blocksize
rng(3);
rows = struct('family', {}, 'method', {}, 'n', {}, 'm', {}, 'L', {}, 'maxrank', {}, ...
              't', {}, 't_repeat', {}, 't_matvec', {}, 'iters', {}, 'res', {}, 'err', {});
fams = {'F1', 'F2', 'F3', 'F4', 'F5'};
for f = 1:numel(fams)
    fam = fams{f};
    for e = 10:opts.maxexp
        n = 2^e;
        switch fam
            case 'F1', G = mp_gen_lrbd(n/2, n, 3, bs);
            case 'F2', G = mp_gen_random(n/2, n, 10, bs);
            case 'F3'
                K = mp_kernel_cauchy(floor(n/3)+1, 3); [L, rb, cb] = mp_tree(K.m, K.n, bs);
                G = mp_hss_kernel(K, L, rb, cb, 1e-9);
            case 'F4'
                K = mp_kernel_conv(n, 2, 'gauss', 2); [L, rb, cb] = mp_tree(K.m, K.n, bs);
                G = mp_hss_kernel(K, L, rb, cb, 1e-13);
            case 'F5'
                m = 3*n/4; xj = sort(((0:m-1)' + 0.5 + 0.5*(rand(m,1)-0.5))/m);
                K = mp_kernel_nudft(xj, n); [L, rb, cb] = mp_tree(K.m, K.n, bs);
                G = mp_hss_kernel(K, L, rb, cb, 1e-10);
        end
        H = mp_hss_from_generators(G);
        [b, xs] = mp_manufactured(H, strcmp(fam, 'F5'));
        [t1, t2, x] = mp_time_solve(H, b, opts.reps);
        tm = mp_timeit(@() H * xs, opts.reps, 1);
        r = newrow(fam, 'ulv', H, G);
        r.t = t1; r.t_repeat = t2; r.t_matvec = tm; r.res = mp_rel(H*x, b); r.err = mp_rel(x, xs);
        rows(end+1) = r; %#ok<AGROW>
        fprintf('%s n=%8d m=%8d L=%2d r=%3d  first H\\b %.3fs  repeat %.4fs  H*x %.4fs | res %.1e err %.1e\n', ...
            fam, r.n, r.m, r.L, r.maxrank, t1, t2, tm, r.res, r.err);
        if strcmp(fam, 'F2') && e <= max(opts.maxdense, opts.maxcg)
            if e <= opts.maxdense                  % dense QR of H' (forms the dense matrix)
                A = full(H);
                [td, ~, xd] = mp_timeit(@() mp_minnorm_dense(A, b), 1, 1);
                r = newrow(fam, 'dense_qr', H, G); r.t = td; r.res = mp_rel(A*xd, b); r.err = mp_rel(xd, xs);
                rows(end+1) = r; %#ok<AGROW>
                clear A
                fprintf('    dense QR %.3fs  err %.1e\n', td, r.err);
            end
            if e <= opts.maxcg                     % CGNE with HSS products, to relative residual 1e-12
                Ht = H';
                t0 = tic; [y, flag, ~, it] = pcg(@(v) H*(Ht*v), b, 1e-12, 5000); xc = Ht*y; tc = toc(t0);
                r = newrow(fam, 'cgne', H, G); r.t = tc; r.iters = it; r.res = mp_rel(H*xc, b); r.err = mp_rel(xc, xs);
                if flag ~= 0, r.method = 'cgne_notconverged'; end
                rows(end+1) = r; %#ok<AGROW>
                fprintf('    CGNE %.3fs (%d its, flag %d)  err %.1e\n', tc, it, flag, r.err);
            end
        end
    end
end
S.rows = rows;
tag = ''; if q, tag = '_quick'; end
mp_write_table(fullfile(here, 'results', ['scaling' tag '.txt']), rows, ...
    {'family', 'method', 'n', 'm', 'L', 'maxrank', 't', 't_repeat', 't_matvec', 'iters', 'res', 'err'}, ...
    sprintf('E3/E4 scaling: t = first H\\b (factor+solve) for ulv; blocksize %d; reps %d', bs, opts.reps));
end

function r = newrow(fam, method, H, G)
mr = 0;
for l = 2:G.L+1
    for i = 1:numel(G.U{l}), mr = max([mr, size(G.U{l}{i},2), size(G.V{l}{i},2)]); end
end
r = struct('family', fam, 'method', method, 'n', size(H, 2), 'm', size(H, 1), 'L', G.L, 'maxrank', mr, ...
           't', NaN, 't_repeat', NaN, 't_matvec', NaN, 'iters', NaN, 'res', NaN, 'err', NaN);
end
