function f = hssp_fro(G)
%HSSP_FRO  Frobenius norm from orthonormalized generators: blocks are disjoint and
% the bases orthonormal, so ||H||_F^2 = sum ||D_i||_F^2 + sum ||B||_F^2.
G = hssp_orth(G);
f2 = sum(cellfun(@(D) norm(D, 'fro')^2, G.D));
for l = 1:G.L
  f2 = f2 + sum(cellfun(@(B) norm(B, 'fro')^2, G.B12{l})) + sum(cellfun(@(B) norm(B, 'fro')^2, G.B21{l}));
end
f = sqrt(f2);
end
