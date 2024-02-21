function H = toeplitz(T)
%TOEPLITZ casts TOEPLITZMAT object to regular matric
    H = toeplitz(T.tc, T.tr);
end

