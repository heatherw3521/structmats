function [L, rb, cb] = cluster_tree(m, n, blocksize)
%CLUSTER_TREE  Row/column cluster boundaries of the tree hss_constructor.m builds.
%   [L, rb, cb] = cluster_tree(m, n, blocksize)
%   L      : leaf level (0 = single dense block)
%   rb{l+1}: 0-based boundaries of the 2^l row clusters at level l (length 2^l+1)
%   cb{l+1}: same for columns.
%   Depth rule copied from hss_constructor.m: the largest depth at which the
%   worst-case (floor) branch still has both dimensions >= blocksize; the
%   first child always gets ceil(size/2).
L = 0; mm = m; nn = n;
while true
    m2 = floor(mm/2); n2 = floor(nn/2);
    if m2 < blocksize || n2 < blocksize
        break
    end
    mm = m2; nn = n2; L = L + 1;
end
rb = cell(L+1,1); cb = cell(L+1,1);
rb{1} = [0 m]; cb{1} = [0 n];
for l = 1:L
    r = 0; c = 0;
    for i = 1:2^(l-1)
        r0 = rb{l}(i); r1 = rb{l}(i+1);
        c0 = cb{l}(i); c1 = cb{l}(i+1);
        r = [r, r0 + ceil((r1-r0)/2), r1]; %#ok<AGROW>
        c = [c, c0 + ceil((c1-c0)/2), c1]; %#ok<AGROW>
    end
    rb{l+1} = r; cb{l+1} = c;
end
end
