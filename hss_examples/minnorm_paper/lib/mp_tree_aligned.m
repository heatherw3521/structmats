function [L, rb, cb] = mp_tree_aligned(rowpos, colpos, blocksize)
%MP_TREE_ALIGNED  Geometry-aligned tree: columns are split exactly as
%   mp_tree splits them; each ROW cut is placed at the first row whose
%   position is >= the position of the corresponding column cut, so that the
%   row and column clusters cover the same piece of the domain.  Positions are
%   real and sorted (e.g. NUDFT sample points and the grid l/N).  Use this when
%   the sampling density is far from uniform (gaps, variable density): the
%   index-based split of hss_constructor.m then misaligns row and column
%   clusters and the off-diagonal ranks grow (used by mp_exp_applications.m).
m = numel(rowpos); n = numel(colpos);
[L, ~, cbt] = mp_tree(m, n, blocksize);
cbL = cbt{L+1};
rbL = zeros(size(cbL));
for i = 2:numel(cbL)-1
    rbL(i) = sum(rowpos(:) < colpos(cbL(i)+1));
end
rbL(end) = m;
rb = cell(L+1,1); cb = cell(L+1,1);
rb{L+1} = rbL; cb{L+1} = cbL;
for l = L-1:-1:0
    rb{l+1} = rb{l+2}(1:2:end); cb{l+1} = cb{l+2}(1:2:end);
end
end
