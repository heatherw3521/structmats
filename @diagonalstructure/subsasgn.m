function h = subsasgn(T, S, b)
%SUBSASGN sets the part of a DIAGONALSTRUCTURE object by referenced subscript.
%   The result may or may not itself be Diagonally structured

% Basic implementation that does not take advantage of structure.
% The way that this will be efficiently overloaded depends on the type
% of diagonally structured matrix
h = subsasgn(full(T), S, b);
end