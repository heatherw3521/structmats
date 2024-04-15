function A = subsref(V,S)
%SUBSREF Subreferencing for VANDERMAT objects.
%   This is a temporary naive implemenation, and only . indexing works.

A = builtin('subsref',V,S);
end

