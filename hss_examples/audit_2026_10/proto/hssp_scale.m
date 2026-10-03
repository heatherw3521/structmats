function G = hssp_scale(G, s)
%HSSP_SCALE  generators of s*H: scale every D and every coupling B (bases unchanged)
for i = 1:numel(G.D), G.D{i} = s*G.D{i}; end
for l = 1:G.L
  for i = 1:numel(G.B12{l}), G.B12{l}{i} = s*G.B12{l}{i}; G.B21{l}{i} = s*G.B21{l}{i}; end
end
end
