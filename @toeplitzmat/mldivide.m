function x = mldivide(T, b)
%MLDIVIDE Left matrix divide for TOEPLITZMAT.
%   Solves the equations Tx = b
x = Toeps(T.tc, T.tr, b);
end

