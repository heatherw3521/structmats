function x = mldivide(H,b)
%MLDIVIDE Left matrix divide for HANKELMAT.
%   Solves Hx = b.
%   We do this efficiently by noting that the flip of a hankelmat
%   is a toeplitzmat! We have a fast solver for toeplitzmat.

error('Hankel Left-divide implementation in progress!')
x = flip(flip(H) \ flip(b));
end

