function x = minnorm(H, b, varargin)
%MINNORM  Minimum-norm or weighted minimum-norm solution of H*x = b.
%   x = minnorm(H, b)               argmin ||x||   s.t. H*x = b  (same as H\b for
%                                   wide H; for square nonsingular H the solution)
%   x = minnorm(H, b, 'Weight', L)  argmin ||L*x|| s.t. H*x = b
%
%   L is either a vector w (L = diag(w), no zero entries) or a cell array
%   {L_1, ..., L_q} of square invertible blocks, one per leaf of H, in leaf
%   order, sized like the leaf column blocks (L = blkdiag(L_1, ..., L_q)).
%   H must have full row rank (wide or square). b may have several columns.
%
%   Method (memo 2): with y = L*x the problem is min ||y|| s.t. (H/L)*y = b,
%   and H/L is HSS on the same tree with only the leaf generators changed
%   (D_tau/L_tau and the leaf column bases). L is applied by solves.
%
%   Factors are kept with H: a repeated call with the same H and the same L
%   (e.g. several right-hand sides, or ADMM/Douglas-Rachford steps) costs one
%   O(n r) solve. A new L refactors and replaces the stored factors. See
%   clearfactors to free the memory.
%
%   See also mldivide, tikhonov, clearfactors.
hss_assertroot(H, 'minnorm');
opts = parse_opts(varargin, {'Weight'}, 'minnorm');
L = opts.Weight;
if H.sz(1) > H.sz(2)
    error('hss:minnorm:tall', ['minnorm needs a wide or square H (got %d x %d); ' ...
        'for least squares use tikhonov.'], H.sz(1), H.sz(2));
end
if isempty(L)
    x = hss_ulvminnormsolve(hss_cachedfactor(H), b);
    return
end
C = H.factorcache;
if ~isempty(C) && ~isempty(C.weighted) && isequal(C.weighted.key, L)
    F = C.weighted.F; cb = C.weighted.cb;
else
    G = hss_to_gen(H);
    cb = G.cb{G.L+1};
    hss_weightops('check', L, cb, 'Weight', true);
    F = hss_ulvminnormsolve(hss_from_gen(hss_gen_weight(G, [], L)));
    if ~isempty(C)
        C.weighted = struct('key', {L}, 'F', F, 'cb', cb);
    end
end
x = hss_weightops('solve', L, hss_ulvminnormsolve(F, b), cb);
end

function opts = parse_opts(args, names, fname)
opts = struct();
for k = 1:numel(names), opts.(names{k}) = []; end
if mod(numel(args), 2) ~= 0
    error(sprintf('hss:%s:options', fname), '%s: options must be name, value pairs.', fname);
end
for k = 1:2:numel(args)
    hit = find(strcmpi(args{k}, names), 1);
    if isempty(hit)
        error(sprintf('hss:%s:options', fname), '%s: unknown option ''%s'' (allowed: %s).', ...
            fname, char(args{k}), strjoin(names, ', '));
    end
    opts.(names{hit}) = args{k+1};
end
end
