function x = mrdivide(T, b)
%MRDIVIDE Right matrix divide for TOEPLITZMAT.
%   Solves the equations xT = b.
%   (note: equivalent to T'x' = (xT)' = b')

%error('Toeplitz right-divide implementation in progress!')
x = (T' \ b')';
end

