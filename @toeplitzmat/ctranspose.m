function v = ctranspose(T)
%CTRANSPOSE for Toeplitz matrices
%   Complex Conjugate Transpose
    v = toeplitzmat(T.tr',T.tc');
end

