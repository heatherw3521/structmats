function G = hssp_orth(G)
%HSSP_ORTH  Make every row basis Ubig and column basis Vbig orthonormal (upward QR
% sweep, R factors pushed into the parent translation and the couplings).
% Also makes the representation proper: ranks <= cluster sizes.
L = G.L;
if L == 0, return; end
for pass = 1:2                       % 1: row bases U, 2: column bases V
  for l = L:-1:1
    R = cell(2^l, 1);
    for i = 1:2^l
      if pass == 1, M = G.U{l+1}{i}; else, M = G.V{l+1}{i}; end
      [Q, R{i}] = qr(M, 0);
      if pass == 1, G.U{l+1}{i} = Q; else, G.V{l+1}{i} = Q; end
    end
    for p = 1:2^(l-1)
      Ra = R{2*p-1}; Rc = R{2*p};
      if pass == 1
        G.B12{l}{p} = Ra * G.B12{l}{p};          % rows of cluster a
        G.B21{l}{p} = Rc * G.B21{l}{p};          % rows of cluster c
        if l >= 2
          T = G.U{l}{p}; ka = size(Ra, 2);
          G.U{l}{p} = [Ra * T(1:ka,:); Rc * T(ka+1:end,:)];
        end
      else
        G.B12{l}{p} = G.B12{l}{p} * Rc';         % columns of cluster c
        G.B21{l}{p} = G.B21{l}{p} * Ra';         % columns of cluster a
        if l >= 2
          T = G.V{l}{p}; ka = size(Ra, 2);
          G.V{l}{p} = [Ra * T(1:ka,:); Rc * T(ka+1:end,:)];
        end
      end
    end
  end
end
end
