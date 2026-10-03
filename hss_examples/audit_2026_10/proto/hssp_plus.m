function G = hssp_plus(G1, G2)
%HSSP_PLUS  generators of H1 + H2 on the same tree (ranks add, no truncation)
if G1.L ~= G2.L || ~isequal(G1.rb, G2.rb) || ~isequal(G1.cb, G2.cb)
  error('hssp_plus:tree', 'H1 and H2 must share the cluster tree');
end
G = G1; L = G.L;
for i = 1:2^L, G.D{i} = G1.D{i} + G2.D{i}; end
if L == 0, return; end
for i = 1:2^L
  G.U{L+1}{i} = [G1.U{L+1}{i}, G2.U{L+1}{i}];
  G.V{L+1}{i} = [G1.V{L+1}{i}, G2.V{L+1}{i}];
end
for l = 1:L-1            % translations of non-leaf nodes at level l
  for i = 1:2^l
    G.U{l+1}{i} = interleave(G1.U{l+1}{i}, G2.U{l+1}{i}, size(G1.U{l+2}{2*i-1},2), size(G2.U{l+2}{2*i-1},2));
    G.V{l+1}{i} = interleave(G1.V{l+1}{i}, G2.V{l+1}{i}, size(G1.V{l+2}{2*i-1},2), size(G2.V{l+2}{2*i-1},2));
  end
end
for l = 0:L-1
  for i = 1:2^l
    G.B12{l+1}{i} = blkdiag(G1.B12{l+1}{i}, G2.B12{l+1}{i});
    G.B21{l+1}{i} = blkdiag(G1.B21{l+1}{i}, G2.B21{l+1}{i});
  end
end
end
function T = interleave(T1, T2, k1a, k2a)
% children bases are [U1_c, U2_c] for c = a, b: rows [a:1, a:2, b:1, b:2]
T = [T1(1:k1a,:), zeros(k1a, size(T2,2));
     zeros(k2a, size(T1,2)), T2(1:k2a,:);
     T1(k1a+1:end,:), zeros(size(T1,1)-k1a, size(T2,2));
     zeros(size(T2,1)-k2a, size(T1,2)), T2(k2a+1:end,:)];
end
