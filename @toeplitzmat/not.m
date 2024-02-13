function h = not(T)
%~  Pointwise logical not for TOEPLITZMAT objects
%   T is a TOEPLITZMAT
%   behaves the same as ~ for matrices in MATLAB

% Check that T is actually logical

if(~islogical(T))
    error( 'TOEPLITZMAT:not:invalidinputtype', ...
        ['Undefined function ''logical not'' for input argument of type %s'], ...
        class(T.tc));
end

% little hack because I'm too lazy to change toepcompare
ftest = @(x) true;
h = toepcompare(T, 1, @(x,y) ~x, "~", ftest);
end