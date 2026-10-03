function [ld, ph] = hssp_logdet(F)
%HSSP_LOGDET  log|det(H)| and phase det(H)/|det(H)| from the stored ULV factors of a
% square H (no size reduction happens for square leaves).
% F is the factor struct stored by H\b: after x = H\b, F = H.factorcache.ulv.
if F.isroot
  if ~strcmp(F.kind, 'lu'), error('hssp_logdet:shape', 'H must be square'); end
  d = diag(F.U); ld = sum(log(abs(d)));
  pv = F.p(:)'; I = eye(numel(pv)); ph = det(I(pv,:)) * prod(d ./ abs(d));
  return
end
ld = 0; ph = 1;
for i = 1:F.nleaves
  if ~isempty(F.Om{i}), error('hssp_logdet:shape', 'H must be square'); end
  d = diag(F.Lt{i}); ld = ld + sum(log(abs(d)));
  ph = ph * prod(d ./ abs(d));
  dP = det(F.P{i}); dQ = det(F.Q{i});     % unit modulus
  ph = ph * (dP/abs(dP)) * (dQ/abs(dQ));
end
[ld2, ph2] = hssp_logdet(F.next);
ld = ld + ld2; ph = ph * ph2;
end
