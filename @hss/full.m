function D = full(H)
% full HSS
D = [blockbuilder(H.A11,H),blockbuilder(H.A12,H);blockbuilder(H.A21,H),blockbuilder(H.A22,H)];
end

