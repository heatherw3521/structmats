function tf = isempty(H)
%ISEMPTY  True if H has no rows or no columns (also for the empty object hss()).
tf = any(size(H) == 0);
end
