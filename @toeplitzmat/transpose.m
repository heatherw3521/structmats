function v = transpose(T)
%TRANSPOSE for TOEPLITZMAT objects
    v = toeplitzmat(T.tr.',T.tc.');
end

