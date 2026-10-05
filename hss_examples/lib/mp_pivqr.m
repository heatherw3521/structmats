function [R, P] = mp_pivqr(A)
%MP_PIVQR  Economy column-pivoted QR, A(:,P) = Q*R, P a permutation vector.
%   Same call in MATLAB (qr(A,'econ','vector')) and GNU Octave (qr(A,0)).
if exist('OCTAVE_VERSION', 'builtin')
    [~, R, P] = qr(A, 0);
else
    [~, R, P] = qr(A, 'econ', 'vector');
end
P = P(:).';
end
