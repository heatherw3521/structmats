function [G, info] = mp_hss_kernel(K, L, rb, cb, tol, opts)
%MP_HSS_KERNEL  Nested interpolative-decomposition HSS construction from entry
%   evaluation, returning GENERATORS (feed to mp_hss_from_generators).
%
%   K : kernel struct (see mp_kernel_cauchy / mp_kernel_conv / mp_kernel_nudft /
%       mp_kernel_toeplitz):  K.m, K.n, K.ent(I,J) (dense block), K.rpos, K.cpos
%       (complex positions, monotone along a curve), K.farr(I,zq), K.farc(zq,J)
%       (bases for the far field through proxy points zq; [] if it vanishes),
%       K.cyclic (unit-circle geometry), K.isreal.
%   L, rb, cb : tree from mp_tree (0-based cluster boundaries).
%   opts.mode = 'proxy' (default) | 'full'
%       'full'  : the rule of hss_constructor.m -- block rows/columns use ALL
%                 other columns/rows (O(m n) entry evaluations).
%       'proxy' : near field explicit (full columns at the leaf level, the
%                 level-(l+1) skeletons of near clusters above it), far field
%                 through opts.nproxy points on a circle of radius
%                 opts.eta_proxy*r around each cluster; near = within
%                 opts.eta_near*r.  O(n) entry evaluations; never forms A.
%   tol : relative ID threshold (|R_jj| >= tol |R_11|).
if nargin < 6, opts = struct(); end
if ~isfield(opts,'mode'), opts.mode = 'proxy'; end
if ~isfield(opts,'eta_near'), opts.eta_near = 3; end
if ~isfield(opts,'eta_proxy'), opts.eta_proxy = 2; end
if ~isfield(opts,'nproxy'), opts.nproxy = 96; end
if ~isfield(opts,'kmax'), opts.kmax = inf; end
m = K.m; n = K.n;
G.L = L; G.rb = rb; G.cb = cb;
nl = 2^L;
G.D = cell(nl,1);
for i = 1:nl
    G.D{i} = K.ent(rb{L+1}(i)+1:rb{L+1}(i+1), cb{L+1}(i)+1:cb{L+1}(i+1));
end
G.U = cell(L+1,1); G.V = cell(L+1,1); G.B12 = cell(L,1); G.B21 = cell(L,1);
info.maxrank = 0;
if L == 0, return, end
zq = exp(2i*pi*((0:opts.nproxy-1)' + 0.5)/opts.nproxy);
srow = cell(L+1,1); scol = cell(L+1,1);
for l = L:-1:1
    nc = 2^l;
    Ul = cell(nc,1); Vl = cell(nc,1); sr = cell(nc,1); sc = cell(nc,1);
    rgeo = clustergeo(K.rpos, rb{l+1}); cgeo = clustergeo(K.cpos, cb{l+1});
    if l < L
        rgeo1 = clustergeo(K.rpos, rb{l+2}); cgeo1 = clustergeo(K.cpos, cb{l+2});
    end
    for i = 1:nc
        if l == L
            cand_r = rb{l+1}(i)+1:rb{l+1}(i+1);
            cand_c = cb{l+1}(i)+1:cb{l+1}(i+1);
        else
            cand_r = [srow{l+2}{2*i-1}, srow{l+2}{2*i}];
            cand_c = [scol{l+2}{2*i-1}, scol{l+2}{2*i}];
        end
        % ---------------- row ID (interactions with columns outside J_i)
        if isempty(cand_r)
            Ul{i} = zeros(0,0); sr{i} = zeros(1,0);
        else
            if strcmp(opts.mode, 'full')
                outside = [1:cb{l+1}(i), cb{l+1}(i+1)+1:n];
                blk = K.ent(cand_r, outside);
            else
                c0 = rgeo.c(i); r0 = rgeo.r(i);
                if isnan(c0), c0 = cgeo.c(i); r0 = cgeo.r(i); end
                if l == L
                    [near, hasfar] = nearclusters(i, nc, cgeo, c0, opts.eta_near*r0, K.cyclic, i);
                    nearcols = idxof(cb{l+1}, near);
                else
                    [near, hasfar] = nearclusters(2*i-1, 2*nc, cgeo1, c0, opts.eta_near*r0, K.cyclic, [2*i-1, 2*i]);
                    nearcols = [scol{l+2}{near}];
                end
                blk = K.ent(cand_r, nearcols);
                if hasfar && ~isempty(K.farr)
                    P = K.farr(cand_r, c0 + opts.eta_proxy*r0*zq);
                    if K.isreal, P = [real(P), imag(P)]; end
                    blk = [blk, P];
                end
            end
            [Z, sel] = mp_id_rows(blk, tol, opts.kmax);
            Ul{i} = Z; sr{i} = cand_r(sel);
        end
        % ---------------- column ID (interactions with rows outside I_i)
        if strcmp(opts.mode, 'full')
            outside = [1:rb{l+1}(i), rb{l+1}(i+1)+1:m];
            blk = K.ent(outside, cand_c);
        else
            c0 = cgeo.c(i); r0 = cgeo.r(i);
            if l == L
                [near, hasfar] = nearclusters(i, nc, rgeo, c0, opts.eta_near*r0, K.cyclic, i);
                nearrows = idxof(rb{l+1}, near);
            else
                [near, hasfar] = nearclusters(2*i-1, 2*nc, rgeo1, c0, opts.eta_near*r0, K.cyclic, [2*i-1, 2*i]);
                nearrows = [srow{l+2}{near}];
            end
            blk = K.ent(nearrows, cand_c);
            if hasfar && ~isempty(K.farc)
                P = K.farc(c0 + opts.eta_proxy*r0*zq, cand_c);
                if K.isreal, P = [real(P); imag(P)]; end
                blk = [blk; P];
            end
        end
        [Y, sel] = mp_id_cols(blk, tol, opts.kmax);
        Vl{i} = Y'; sc{i} = cand_c(sel);
        info.maxrank = max([info.maxrank, size(Ul{i},2), size(Vl{i},2)]);
    end
    G.U{l+1} = Ul; G.V{l+1} = Vl; srow{l+1} = sr; scol{l+1} = sc;
end
for l = 0:L-1
    G.B12{l+1} = cell(2^l,1); G.B21{l+1} = cell(2^l,1);
    for i = 1:2^l
        a = 2*i-1; c = 2*i;
        G.B12{l+1}{i} = K.ent(srow{l+2}{a}, scol{l+2}{c});
        G.B21{l+1}{i} = K.ent(srow{l+2}{c}, scol{l+2}{a});
    end
end
info.srow = srow; info.scol = scol;
end

function g = clustergeo(pos, bnd)
nc = numel(bnd) - 1;
g.c = nan(nc,1); g.r = zeros(nc,1);
for i = 1:nc
    p = pos(bnd(i)+1:bnd(i+1));
    if isempty(p), continue, end
    c = complex((min(real(p))+max(real(p)))/2, (min(imag(p))+max(imag(p)))/2);
    g.c(i) = c; g.r(i) = max(max(abs(p - c)), 1e-300);
end
end

function [near, hasfar] = nearclusters(i, nclus, geo, c, rad, cyclic, exclude)
% clusters (not in exclude) whose bounding disks meet |z - c| < rad, scanning
% outward from cluster i along the 1-D ordering (1-based indices)
isnear = false(nclus,1);
for dirn = [1, -1]
    for step = 1:nclus-1
        j = i + dirn*step;
        if cyclic
            j = mod(j-1, nclus) + 1;
        elseif j < 1 || j > nclus
            break
        end
        if any(j == exclude), continue, end
        if isnan(geo.c(j)), isnear(j) = true; continue, end
        if abs(geo.c(j) - c) - geo.r(j) < rad
            isnear(j) = true;
        else
            break
        end
    end
end
near = find(isnear).';
hasfar = (nclus - numel(near) - numel(exclude)) > 0;
end

function idx = idxof(bnd, clusters)
idx = zeros(1,0);
for j = clusters
    idx = [idx, bnd(j)+1:bnd(j+1)]; %#ok<AGROW>
end
end
