function B = hssp_ulv_adjoint(F, Y)
%HSSP_ULV_ADJOINT  Adjoint of the stored min-norm solve: if x = S(b) is the solve with the
% factors F of a wide/square full-row-rank W (S = pinv(W)), this returns
% S' * Y = pinv(W)' * Y = pinv(W') * Y.  With W = H' for a tall H of full
% column rank, this is the least-squares solution of H x = Y.
% F is the factor struct stored by H\b: after x = H\b, F = H.factorcache.ulv.
if F.isroot
  if strcmp(F.kind, 'lu')          % X = U \ (L \ B(p,:))
    Z = F.L' \ (F.U' \ Y);
    B = zeros(size(Z), 'like', Z); B(F.p, :) = Z;
  else                             % X = V * ((U'*B) ./ s)
    B = F.U * ((F.V' * Y) ./ F.s);
  end
  return
end
ns = size(Y, 2);
% adjoint of the reconstruction X_i = Om_i' Q_i' W_i
Wt = zeros(F.co(end), ns);
for i = 1:F.nleaves
  Yi = Y(F.xo(i) + (1:F.n(i)), :);
  if ~isempty(F.Om{i}), Yi = F.Om{i} * Yi; end
  Wt(F.co(i) + (1:F.r(i)), :) = F.Q{i} * Yi;
end
% adjoint of the next level, scattered back to the kept rows
g = zeros(F.ro(end), ns);
g(F.keep_rows, :) = hssp_ulv_adjoint(F.next, Wt(F.free_cols, :));
% a = z_forced - S_f' * HULV' * g
a = Wt - F.HULV' * g;              % only the forced entries are used below
B1 = g;
for i = 1:F.nleaves
  t = F.t(i);
  if t > 0
    B1(F.ro(i) + (1:t), :) = B1(F.ro(i) + (1:t), :) + F.Lt{i}' \ a(F.co(i) + (1:t), :);
  end
end
% adjoint of the row compression b_i <- P_i' b_i
B = zeros(size(B1), 'like', B1);
for i = 1:F.nleaves
  rows = F.ro(i) + (1:F.m(i));
  B(rows, :) = F.P{i} * B1(rows, :);
end
end
