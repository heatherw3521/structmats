function h = subsasgn(T, S, b)
%SUBSASGN sets the part of a Toeplitzmat object by referenced subscript.
%   The result may or may not itself be Toeplitz

    h = subsasgn(toeplitz(T), S, b);

%% TO DO: Make this function smart enough to detect when the result
%   will be Toeplitz. Not sure how useful this will be.

end