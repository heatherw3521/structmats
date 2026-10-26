function S = extract(H, I, J)
%EXTRACT  S = H(I, J) for global index vectors I, J (any order, repeats
%   allowed) without forming any off-diagonal block. Called by subsref after
%   the indices are checked.
%
%   Every off-diagonal block is Ubig * B * Ybig, with nested bases: the row
%   basis of a cluster is blkdiag(row bases of its two children) times a
%   translation matrix, down to the leaf bases. Only the rows of the row
%   bases at I and the columns of the column bases at J are computed, bottom
%   up through the tree, so the cost is O((|I| + |J|) r log n + |I| |J| r)
%   instead of O(n^2) for a single entry.
I = I(:); J = J(:);
if H.isleaf                  % a single-leaf root: just the dense block
    S = H.D(I, J);
    return
end
S = node(H, I, J);
end


function [S, Lr, Rc] = node(D, I, J)
% D: diagonal node; I, J: global indices inside D.
% S  = D(I, J)
% Lr = rows I of blkdiag(Ubig of D's two children)      (|I| x (ka + kc))
% Rc = columns J of blkdiag(Ybig of D's two children)   ((ka + kc) x |J|)
% (Lr, Rc are used by D's parent to form D's own bases; empty for a leaf.)
if D.isleaf
    S = D.D(I - D.Ir(1) + 1, J - D.Ic(1) + 1);
    Lr = []; Rc = [];
    return
end
if isempty(I) && isempty(J)              % nothing requested below D
    S = zeros(0, 0);
    Lr = zeros(0, size(D.A12.Z, 2) + size(D.A21.Z, 2));
    Rc = zeros(size(D.A21.Y, 1) + size(D.A12.Y, 1), 0);
    return
end
a = D.A11; c = D.A22;
ia = I <= a.Ir(2); ic = ~ia;              % rows in the first / second child
ja = J <= a.Ic(2); jc = ~ja;              % columns likewise

[Sa, Lra, Rca] = node(a, I(ia), J(ja));
[Sc, Lrc, Rcc] = node(c, I(ic), J(jc));

% bases of the two children restricted to the requested rows / columns:
%   rows of a use D.A12.Z (translation, or the leaf basis at the leaf level),
%   rows of c use D.A21.Z; columns of a use D.A21.Y, columns of c D.A12.Y
Ua = rowbasis(a, I(ia), Lra, D.A12.Z);
Uc = rowbasis(c, I(ic), Lrc, D.A21.Z);
Ya = colbasis(a, J(ja), Rca, D.A21.Y);
Yc = colbasis(c, J(jc), Rcc, D.A12.Y);

S = zeros(numel(I), numel(J), 'like', Sa);
S(ia, ja) = Sa;
S(ic, jc) = Sc;
if any(ia) && any(jc)
    S(ia, jc) = Ua * D.A12.lrcomponent * Yc;
end
if any(ic) && any(ja)
    S(ic, ja) = Uc * D.A21.lrcomponent * Ya;
end

if nargout > 1
    ka = size(D.A12.Z, 2); kc = size(D.A21.Z, 2);
    Lr = zeros(numel(I), ka + kc, 'like', S);
    Lr(ia, 1:ka) = Ua;
    Lr(ic, ka+1:end) = Uc;
    la = size(D.A21.Y, 1); lc = size(D.A12.Y, 1);
    Rc = zeros(la + lc, numel(J), 'like', S);
    Rc(1:la, ja) = Ya;
    Rc(la+1:end, jc) = Yc;
end
end


function U = rowbasis(child, I, Lr, Z)
% rows I of the child's row basis: the leaf basis Z itself at a leaf,
% otherwise blkdiag(grandchildren bases)(I,:) * translation Z
if isempty(I)
    U = zeros(0, size(Z, 2), 'like', Z);
elseif child.isleaf
    U = Z(I - child.Ir(1) + 1, :);
else
    U = Lr * Z;
end
end


function Y = colbasis(child, J, Rc, Ym)
% columns J of the child's column basis (stored as Y = V')
if isempty(J)
    Y = zeros(size(Ym, 1), 0, 'like', Ym);
elseif child.isleaf
    Y = Ym(:, J - child.Ic(1) + 1);
else
    Y = Ym * Rc;
end
end
