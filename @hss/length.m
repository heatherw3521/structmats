function L = length(H)
%LENGTH  Largest dimension of H (0 if H is empty), as for a dense matrix.
s = size(H);
if any(s == 0), L = 0; else, L = max(s); end
end
