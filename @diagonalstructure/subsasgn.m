function h = subsasgn(T, S, b)
%SUBSASGN sets the part of a DIAGONALSTRUCTURE object by referenced subscript.
%   The result may or may not itself be Toeplitz

    h = subsasgn(full(T), S, b);

%% TO DO: Make this function smart enough to detect when the result
%   will have the same structure

end