function e = end(H, k, n)
%END  Last index in H(..., end, ...): size(H, k), or numel for linear indexing.
s = H.sz;
if n == 1
    e = prod(s);
elseif k <= 2
    e = s(k);
else
    e = 1;
end
end
