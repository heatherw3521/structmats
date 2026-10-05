function G = gen_lowrank_blockdiag(m, n, k, blocksize, seed)
%GEN_LOWRANK_BLOCKDIAG  Family F1: low rank plus block diagonal, in generator form.
%   A = X*Y + blkdiag(E_tau),  X (m x k), Y (k x n), E_tau (leaf diagonal blocks),
%   all entries iid U[0,1].  Off-diagonal blocks have rank exactly k at every
%   level: U_tau = X(I_tau,:), V_tau = Y(:,J_tau).', R = W = [I;I], B = I.
%   Never forms the dense matrix (O(n k) memory).  Tree = cluster_tree(m,n,blocksize).
if nargin >= 5, rng(seed); end
[L, rb, cb] = cluster_tree(m, n, blocksize);
X = rand(m, k); Y = rand(k, n);
G.L = L; G.rb = rb; G.cb = cb; G.blocksize = blocksize;
G.X = X; G.Y = Y;   % kept for verification only (O(n k))
nl = 2^L;
G.D = cell(nl, 1);
for i = 1:nl
    r = rb{L+1}(i)+1:rb{L+1}(i+1); c = cb{L+1}(i)+1:cb{L+1}(i+1);
    G.D{i} = X(r,:)*Y(:,c) + rand(numel(r), numel(c));
end
G.U = cell(L+1,1); G.V = cell(L+1,1); G.B12 = cell(L,1); G.B21 = cell(L,1);
if L == 0, return, end
G.U{L+1} = cell(nl,1); G.V{L+1} = cell(nl,1);
for i = 1:nl
    G.U{L+1}{i} = X(rb{L+1}(i)+1:rb{L+1}(i+1), :);
    G.V{L+1}{i} = Y(:, cb{L+1}(i)+1:cb{L+1}(i+1)).';
end
I = eye(k);
for l = 1:L-1
    G.U{l+1} = repmat({[I; I]}, 2^l, 1);
    G.V{l+1} = repmat({[I; I]}, 2^l, 1);
end
for l = 0:L-1
    G.B12{l+1} = repmat({I}, 2^l, 1);
    G.B21{l+1} = repmat({I}, 2^l, 1);
end
end
