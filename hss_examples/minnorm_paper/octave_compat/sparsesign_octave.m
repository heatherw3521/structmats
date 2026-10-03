function S = sparsesign(d, m, zeta)
% Octave test-harness replacement for the sparsesign MEX (same distribution:
% d x m sparse sign embedding, zeta nonzeros per column, values +-1/sqrt(zeta)).
zeta = min(zeta, d);
rows = zeros(m*zeta,1);
for i = 1:m
  rows((i-1)*zeta+1:i*zeta) = randperm(d, zeta);
end
cols = kron((1:m)', ones(zeta,1));
vals = (2*randi(2, m*zeta, 1) - 3) / sqrt(zeta);
S = sparse(rows, cols, vals, d, m);
end
