function A = subsasgn(V,S,B)
%SUBSASGN Subassign for VANDERMAT objects.
%   Naive temporary implementation below.

A = builtin('subsref',V,S,B);
end

