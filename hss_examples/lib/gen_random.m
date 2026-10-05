function G = gen_random(m, n, k, blocksize, seed, opts)
%GEN_RANDOM  Family F2: exact random HSS with rank k on every node.
%   D_tau iid N(0,1)/sqrt(p_tau); U, V, R, W with orthonormal columns (QR of
%   Gaussians); B iid N(0,1)*coupling/sqrt(k).  opts.complex, opts.coupling.
%   Never forms the dense matrix.
if nargin < 6, opts = struct(); end
if ~isfield(opts, 'complex'), opts.complex = false; end
if ~isfield(opts, 'coupling'), opts.coupling = 1; end
if nargin >= 5 && ~isempty(seed), rng(seed); end
g = @(a,b) randn(a,b) + 1i*opts.complex*randn(a,b);
[L, rb, cb] = cluster_tree(m, n, blocksize);
G.L = L; G.rb = rb; G.cb = cb; G.blocksize = blocksize;
nl = 2^L;
G.D = cell(nl,1);
for i = 1:nl
    nr = rb{L+1}(i+1)-rb{L+1}(i); nc = cb{L+1}(i+1)-cb{L+1}(i);
    G.D{i} = g(nr, nc) / sqrt(nc);
end
G.U = cell(L+1,1); G.V = cell(L+1,1); G.B12 = cell(L,1); G.B21 = cell(L,1);
if L == 0, return, end
for l = 1:L
    G.U{l+1} = cell(2^l,1); G.V{l+1} = cell(2^l,1);
    for i = 1:2^l
        if l == L
            a = rb{L+1}(i+1)-rb{L+1}(i); b = cb{L+1}(i+1)-cb{L+1}(i);
        else
            a = 2*k; b = 2*k;
        end
        [Q,~] = qr(g(a,k), 0); G.U{l+1}{i} = Q;
        [Q,~] = qr(g(b,k), 0); G.V{l+1}{i} = Q;
    end
end
for l = 0:L-1
    G.B12{l+1} = cell(2^l,1); G.B21{l+1} = cell(2^l,1);
    for i = 1:2^l
        G.B12{l+1}{i} = g(k,k)*opts.coupling/sqrt(k);
        G.B21{l+1}{i} = g(k,k)*opts.coupling/sqrt(k);
    end
end
end
