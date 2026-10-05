function e = mp_rel(a, b)
%MP_REL  relative error ||a - b|| / ||b||
e = norm(a - b) / norm(b);
end
