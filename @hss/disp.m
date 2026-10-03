function disp(H)
%DISP  One-line summary of an HSS matrix (never the whole tree).
if isempty(H.sz)
    fprintf('  empty hss object\n\n');
    return
end
if isempty(H.isroot) || ~H.isroot
    fprintf('  piece of an HSS matrix: %dx%d block, rows %d-%d, columns %d-%d of the full matrix\n\n', ...
        H.sz(1), H.sz(2), H.Ir(1), H.Ir(2), H.Ic(1), H.Ic(2));
    return
end
[kmax, nleaf] = walk(H, 0, 0);
cached = ~isempty(H.factorcache) && ~isempty(H.factorcache.ulv);
fprintf('  %dx%d HSS matrix: %d levels, %d leaves, max off-diagonal rank %d', ...
    H.sz(1), H.sz(2), H.levelcount, nleaf, kmax);
if cached, fprintf(', ULV factors stored'); end
fprintf('\n\n');
end

function [kmax, nleaf] = walk(N, kmax, nleaf)
if N.isleaf
    nleaf = nleaf + 1;
    return
end
kmax = max([kmax, size(N.A12.lrcomponent), size(N.A21.lrcomponent)]);
[kmax, nleaf] = walk(N.A11, kmax, nleaf);
[kmax, nleaf] = walk(N.A22, kmax, nleaf);
end
