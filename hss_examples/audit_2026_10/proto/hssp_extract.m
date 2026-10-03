function S = hssp_extract(G, I, J)
%HSSP_EXTRACT  S = H(I, J) from the generators without building any off-diagonal
%   block: cost O((|I| + |J|) r log n + |I||J|) instead of O(n^2) per call.
% I, J: integer vectors (any order, repeats allowed), validated.
m = G.rb{1}(end); n = G.cb{1}(end); L = G.L;
I = I(:); J = J(:);
if any(I < 1 | I > m | I ~= round(I)) || any(J < 1 | J > n | J ~= round(J))
  error('hssp_extract:index', 'index out of range or not an integer');
end
S = zeros(numel(I), numel(J));
% leaf of every requested row / column
li = 1 + sum(bsxfun(@gt, I, G.rb{L+1}(2:end-1)), 2);
lj = 1 + sum(bsxfun(@gt, J, G.cb{L+1}(2:end-1)), 2);
nl = 2^L;
Ip = cell(nl,1); Jp = cell(nl,1);   % positions into I/J, per leaf
for i = 1:nl, Ip{i} = find(li == i); Jp{i} = find(lj == i); end
% diagonal blocks
for i = 1:nl
  if ~isempty(Ip{i}) && ~isempty(Jp{i})
    S(Ip{i}, Jp{i}) = G.D{i}(I(Ip{i}) - G.rb{L+1}(i), J(Jp{i}) - G.cb{L+1}(i));
  end
end
if L == 0, return; end
% partial bases restricted to the requested rows/columns, bottom up
Up = cell(nl,1); Vp = cell(nl,1); Ipos = Ip; Jpos = Jp;
for i = 1:nl
  Up{i} = G.U{L+1}{i}(I(Ip{i}) - G.rb{L+1}(i), :);
  Vp{i} = G.V{L+1}{i}(J(Jp{i}) - G.cb{L+1}(i), :);
end
for l = L-1:-1:0
  for p = 1:2^l
    a = 2*p-1; c = 2*p;
    % couplings between the children of (l, p)
    if ~isempty(Ipos{a}) && ~isempty(Jpos{c})
      S(Ipos{a}, Jpos{c}) = Up{a} * G.B12{l+1}{p} * Vp{c}';
    end
    if ~isempty(Ipos{c}) && ~isempty(Jpos{a})
      S(Ipos{c}, Jpos{a}) = Up{c} * G.B21{l+1}{p} * Vp{a}';
    end
  end
  if l == 0, break; end
  Un = cell(2^l,1); Vn = Un; In = Un; Jn = Un;
  for p = 1:2^l
    a = 2*p-1; c = 2*p;
    Un{p} = blkdiag(Up{a}, Up{c}) * G.U{l+1}{p};
    Vn{p} = blkdiag(Vp{a}, Vp{c}) * G.V{l+1}{p};
    In{p} = [Ipos{a}; Ipos{c}]; Jn{p} = [Jpos{a}; Jpos{c}];
  end
  Up = Un; Vp = Vn; Ipos = In; Jpos = Jn;
end
end
