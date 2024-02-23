function h = subsref(T,S)
%SUBSREF gets the part of a Toeplitzmat object by referenced subscript.
%   The result may or may not itself be Toeplitz

    h = subsref(toeplitz(T), S);

%% TO DO: Make this function smart enough to detect when the result
%   will be Toeplitz

end

