function H = uminus(H)
%UMINUS  -H for an HSS matrix.
hss_assertroot(H, '-H');
H = hss_scale(H, -1);
end
