function Y = hssp_matvec(G, X)
%HSSP_MATVEC  Y = H*X from the generators with per-level cell arrays
%   (upward pass for the column coefficients, downward pass for the row
%   coefficients). Same arithmetic as @hss/private/hss_matvec.m, but no
%   dictionary: plain cell arrays indexed by (level, node).
L = G.L; ns = size(X, 2);
if size(X, 1) ~= G.cb{1}(end)
  error('hssp_matvec:dimension', 'X has %d rows; H has %d columns.', size(X,1), G.cb{1}(end));
end
if L == 0, Y = G.D{1} * X; return; end
cb = G.cb{L+1}; rb = G.rb{L+1};
xc = cell(L+1, 1); xc{L+1} = cell(2^L, 1);
for i = 1:2^L, xc{L+1}{i} = G.V{L+1}{i}' * X(cb(i)+1:cb(i+1), :); end
for l = L-1:-1:1
  xc{l+1} = cell(2^l, 1);
  for p = 1:2^l, xc{l+1}{p} = G.V{l+1}{p}' * [xc{l+2}{2*p-1}; xc{l+2}{2*p}]; end
end
yc = cell(L+1, 1); yc{1} = {[]};
for l = 0:L-1
  yc{l+2} = cell(2^(l+1), 1);
  for p = 1:2^l
    a = 2*p-1; c = 2*p;
    ya = G.B12{l+1}{p} * xc{l+2}{c};  yb = G.B21{l+1}{p} * xc{l+2}{a};
    if l >= 1
      T = G.U{l+1}{p}; ka = size(G.U{l+2}{a}, 2); z = T * yc{l+1}{p};
      ya = ya + z(1:ka, :); yb = yb + z(ka+1:end, :);
    end
    yc{l+2}{a} = ya; yc{l+2}{c} = yb;
  end
end
Y = zeros(rb(end), ns, 'like', X(1) * G.D{1}(1));
for i = 1:2^L
  Y(rb(i)+1:rb(i+1), :) = G.D{i} * X(cb(i)+1:cb(i+1), :) + G.U{L+1}{i} * yc{L+1}{i};
end
end
