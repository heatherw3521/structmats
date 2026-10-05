function G = hssp_compress(G, tol)
%HSSP_COMPRESS  Recompress HSS generators: orthonormalize, then truncate every row and
% column basis top-down by an SVD of its full off-diagonal coefficient.
% Keeps singular values > tol * (largest singular value of that block row).
if nargin < 2, tol = 1e-12; end
G = hssp_orth(G);
L = G.L;
if L == 0, return; end
for pass = 1:2
  C = {[]};                          % compact coefficient of the root (none)
  for l = 1:L
    Cn = cell(2^l, 1);
    for i = 1:2^l
      p = ceil(i/2); first = mod(i,2) == 1;
      % coupling of cluster i with its sibling, seen from i's side
      if pass == 1
        if first, Bs = G.B12{l}{p}; else, Bs = G.B21{l}{p}; end
      else
        if first, Bs = G.B21{l}{p}'; else, Bs = G.B12{l}{p}'; end
      end
      S = Bs;
      if l >= 2
        if pass == 1, T = G.U{l}{p}; else, T = G.V{l}{p}; end
        ka = size_child(G, pass, l, 2*p-1);
        if first, Tb = T(1:ka,:); else, Tb = T(ka+1:end,:); end
        S = [S, Tb * C{p}];
      end
      if isempty(S), k = 0; W = zeros(size(S,1),0);
      else
        [W, Sig, ~] = svd(S, 'econ'); s = diag(Sig);
        if isempty(s) || s(1) == 0, k = 0; else, k = sum(s > tol*s(1)); end
        W = W(:, 1:k);
      end
      Cn{i} = W' * S;
      % apply W: basis/translation of i, coupling with sibling, parent translation rows
      if pass == 1
        G.U{l+1}{i} = G.U{l+1}{i} * W;
        if first, G.B12{l}{p} = W' * G.B12{l}{p}; else, G.B21{l}{p} = W' * G.B21{l}{p}; end
      else
        G.V{l+1}{i} = G.V{l+1}{i} * W;
        if first, G.B21{l}{p} = G.B21{l}{p} * W; else, G.B12{l}{p} = G.B12{l}{p} * W; end
      end
      if l >= 2
        if pass == 1, T = G.U{l}{p}; else, T = G.V{l}{p}; end
        ka = size_child_rows(T, G, pass, l, p, first, W);
        if first
          T = [W' * T(1:ka,:); T(ka+1:end,:)];
        else
          T = [T(1:ka,:); W' * T(ka+1:end,:)];
        end
        if pass == 1, G.U{l}{p} = T; else, G.V{l}{p} = T; end
      end
    end
    C = Cn;
  end
end
end

function k = size_child(G, pass, l, j)
% current rank (number of columns) of the basis of node (l, j)
if pass == 1, k = size(G.U{l+1}{j}, 2); else, k = size(G.V{l+1}{j}, 2); end
end

function ka = size_child_rows(T, G, pass, l, p, first, W)
% number of rows of parent translation T belonging to the first child
% (the first child may already have been truncated in this sweep)
if first
  ka = size(W, 1);                      % old rank of the first child
else
  ka = size_child(G, pass, l, 2*p-1);   % first child already truncated
end
end
