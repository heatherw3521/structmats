function out = hss_weightops(op, W, varargin)
%HSS_WEIGHTOPS  Checks and products for leaf-conforming weights.
%   W is [] (identity), a vector (diagonal) or a cell array of square blocks,
%   one per leaf in leaf order; bnd are the 0-based leaf boundaries
%   (G.rb{end} for row weights, G.cb{end} for column weights).
%     hss_weightops('check', W, bnd, name, invertible)   errors unless W fits
%     hss_weightops('apply', W, X, bnd)      W * X
%     hss_weightops('solve', W, X, bnd)      W \ X (blockwise solves)
switch op
    case 'check'
        bnd = varargin{1}; name = varargin{2}; inv_req = varargin{3}; n = bnd(end);
        if isempty(W), out = true; return; end
        if iscell(W)
            if numel(W) ~= numel(bnd) - 1
                error('hss:weight:blocks', '%s has %d blocks; H has %d leaves.', ...
                    name, numel(W), numel(bnd) - 1);
            end
            for i = 1:numel(W)
                ni = bnd(i+1) - bnd(i);
                if ~isequal(size(W{i}), [ni ni])
                    error('hss:weight:blocks', '%s{%d} is %dx%d; leaf %d has %d %s.', ...
                        name, i, size(W{i}, 1), size(W{i}, 2), i, ni, 'entries');
                end
                if ~all(isfinite(W{i}(:)))
                    error('hss:weight:nonFinite', '%s{%d} contains NaN or Inf.', name, i);
                end
                % a singular block would make the weighted problem meaningless
                if inv_req && ni > 0 && ~(rcond(full(W{i})) > eps)
                    error('hss:weight:singular', ['%s{%d} is singular to working precision ' ...
                        '(rcond = %.1e); it must be invertible.'], name, i, rcond(full(W{i})));
                end
            end
        else
            if ~isvector(W) || numel(W) ~= n
                error('hss:weight:size', '%s must be [], a vector of length %d, or a cell of leaf blocks.', name, n);
            end
            if ~all(isfinite(W(:)))
                error('hss:weight:nonFinite', '%s contains NaN or Inf.', name);
            end
            if inv_req && any(W(:) == 0)
                error('hss:weight:singular', '%s has zero entries; it must be invertible.', name);
            end
        end
        out = true;
    case {'apply', 'solve'}
        X = varargin{1}; bnd = varargin{2};
        if isempty(W), out = X; return; end
        if ~iscell(W)
            if strcmp(op, 'apply'), out = W(:) .* X; else, out = X ./ W(:); end
            return
        end
        out = zeros(size(X), 'like', X);      % becomes complex on assignment if needed
        for i = 1:numel(W)
            r = bnd(i)+1:bnd(i+1);
            if strcmp(op, 'apply'), out(r, :) = W{i} * X(r, :); else, out(r, :) = W{i} \ X(r, :); end
        end
    otherwise
        error('hss_weightops: unknown op %s', op);
end
end
