function v = ctranspose(T)
%CTRANSPOSE Conjugate Transpose for TOEPLITZMAT objects
    v = toeplitzmat(T.tr',T.tc');
end

