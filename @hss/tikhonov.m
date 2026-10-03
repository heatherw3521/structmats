function [x, r] = tikhonov(H, b, lambda, varargin)
%TIKHONOV  Tikhonov-regularized solve with an HSS matrix of any shape.
%   x = tikhonov(H, b, lambda)
%       argmin ||H*x - b||^2 + lambda^2 ||x||^2
%   x = tikhonov(H, b, lambda, 'Weight', L, 'DataWeight', S)
%       argmin ||S*(H*x - b)||^2 + lambda^2 ||L*x||^2
%   [x, r] = tikhonov(...) also returns the weighted residual r = S*(H*x - b).
%
%   lambda > 0. H may be wide, square or tall. L and S are each a vector
%   (diagonal) or a cell array of square leaf blocks in leaf order: L's blocks
%   sized like H's leaf column blocks (L invertible), S's like its leaf row
%   blocks. Omitted: identity. b may have several columns.
%
%   Method (memo 2): standard form Ht = S*H/L, bt = S*b; then the minimum-norm
%   solution of the augmented wide system [Ht, lambda*I] z = bt, with the
%   identity columns interleaved leaf by leaf, is HSS on the same tree and has
%   full row rank for every lambda > 0. x = L \ z_x, r = -lambda*z_s.
%
%   Factors are kept with H for the last (lambda, L, S): repeating the call
%   with the same values (new right-hand sides) costs one O(n r) solve. Each
%   new lambda needs a new factorization (a lambda sweep of length K costs K).
%
%   See also minnorm, mldivide, clearfactors.
hss_assertroot(H, 'tikhonov');
opts = parse_opts(varargin, {'Weight', 'DataWeight'}, 'tikhonov');
L = opts.Weight; S = opts.DataWeight;
if ~(isscalar(lambda) && isreal(lambda) && lambda > 0)
    error('hss:tikhonov:lambda', 'lambda must be a real scalar > 0.');
end
key = {lambda, L, S};
C = H.factorcache;
if ~isempty(C) && ~isempty(C.tikhonov) && isequal(C.tikhonov.key, key)
    T = C.tikhonov;
else
    G = hss_to_gen(H);
    T.rb = G.rb{G.L+1}; T.cb = G.cb{G.L+1};
    hss_weightops('check', L, T.cb, 'Weight', true);
    hss_weightops('check', S, T.rb, 'DataWeight', false);
    [Ga, T.ix, T.is] = hss_gen_augment(hss_gen_weight(G, S, L), lambda);
    T.F = hss_ulvminnormsolve(hss_from_gen(Ga));
    T.key = key;
    if ~isempty(C)
        C.tikhonov = T;
    end
end
z = hss_ulvminnormsolve(T.F, hss_weightops('apply', S, b, T.rb));
x = hss_weightops('solve', L, z(T.ix, :), T.cb);
if nargout > 1
    r = -lambda * z(T.is, :);
end
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
