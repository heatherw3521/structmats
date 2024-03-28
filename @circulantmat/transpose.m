function C = transpose(C)
%TRANSPOSE Transpose operation for CIRCULANTMAT objects

    tc = C.tc;
    C = circulantmat([tc(1); ...
                     flip(tc(2:end))]);
end

