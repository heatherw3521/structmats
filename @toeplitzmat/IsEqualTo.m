function h = IsEqualTo(T, H)
%ISEQUALTO Strict Equals function for TOEPLITZMAT

if(isa(T, 'toeplitzmat') && isa(H, 'toeplitzmat'))
    h = (isequal(T.tc,H.tc) && isequal(T.tr,H.tr));
elseif(isa(H, 'toeplitzmat'))
    h = isequal(H,T);
else
    h = isequal(full(T),full(H));
end

end

