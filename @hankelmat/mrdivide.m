function x = mrdivide(H, b)
%MRDIVIDE Right matrix divide for HANKELMAT.
%   Solves the equations xH = b.
%   (note: equivalent to H'x'=(xH)' = b')
x = (H' \ b')';
end

