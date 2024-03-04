function res_toeps = structure_retained(T, I, J)
%STRUCTURE_RETAINED determines whether the indexing I,J corresponds to a
%   toeplitz submatrix of a toeplitz matrix.
%
%   The idea here is that the result will be toeplitz if the difference
%   matrix I-J is toeplitz. We compute whether the structure is retained
%   without forming this matrix, as doing so would take up to O(nm) space.
%
%   In fact, the problem reduces the whether or not I,J are uniformly
%   spaced with same spacing! For instance, two examples of
%   structure-retaining indexings are:
%
%   1. I = [1 2 3 4]', J = [5 6 7]
%   2. I = [1 4 7 10]', J = [3 6 9 12 15 18]
%
%   The case of wrap-arounds is slightly more involved, and has yet to be
%   implemented. This is not a priority, as the code will still work for
%   these inputs, it just won't detect that the result should be toeplitz!

% Helper function to determine whether an indexing is spaced uniformly
spaced = @(k, dk) isequal(k, k(1) + dk*(0:length(k)-1));

% Compute whether the result is toeplitz
if(isequal(I,':') && isequal(J,':'))
    res_toeps = true;
elseif(isequal(I,':')) % Is J consecutive
    res_toeps = spaced(J, 1);
elseif(isequal(J,':')) % Is I consecutive?
    res_toeps = spaced(I,1);
elseif(length(I) == 1 || length(J) == 1) % Result is a row or a vector 
    res_toeps = true;
else
    di = I(2)-I(1);
    res_toeps = spaced(J,di);
end

end

