function v = transpose(T)
%TRANSPOSE for Toeplitz matrices
    v = toeplitzmat(T.tr.',T.tc.');
end

