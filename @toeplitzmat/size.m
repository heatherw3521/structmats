function s = size(T)
%SIZE for TOEPLITZMAT objects
%  returns the size of toeplitz matrix
    s = [size(T.tc, 1), size(T.tr, 2)];
end