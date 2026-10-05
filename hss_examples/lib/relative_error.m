function e = relative_error(a, b)
%RELATIVE_ERROR  relative error ||a - b|| / ||b||
e = norm(a - b) / norm(b);
end
