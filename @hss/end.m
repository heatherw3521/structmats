function e = end(H, k, n)
%END  Last index in H(..., end, ...): size(H, k), or numel for linear indexing.
%   Without this method MATLAB used the size of the 1x1 object, so end was
%   always 1 and H(end,end) returned H(1,1) (audit B02).
s = H.sz;
if n == 1
    e = prod(s);
elseif k <= 2
    e = s(k);
else
    e = 1;
end
end
