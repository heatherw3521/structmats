function H = toeplitz(T)
%TOEPLITZ casts TOEPLITZMAT object to regular matrix

    if(islogical(T.tc))
        H = logical(toeplitz(double(T.tc), double(T.tr)));
    else
        H = toeplitz(T.tc, T.tr);
    end
end

