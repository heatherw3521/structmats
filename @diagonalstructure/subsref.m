function h = subsref(T, S)
%SUBSREF gets the part of a DIAGONALSTRUCTURE object by referenced subscript.
%   The result may or may not itself be Diagonally structured

% Basic implementation that does not take advantage of structure.
% The way that this will be efficiently overloaded depends on the type
% of diagonally structured matrix
h = subsref(full(T), S);

end