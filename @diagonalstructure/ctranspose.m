function v = ctranspose(T)
%CTRANSPOSE Conjugate Transpose for DIAGONALSTRUCTURE objects

    % Quick hack to make the most of diag_apply!
    ftest = @(x) true;
    v = diag_apply(@(x,y) conj(x), ftest, "conj", T.', 1);
end

