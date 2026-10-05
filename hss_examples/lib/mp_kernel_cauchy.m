function K = mp_kernel_cauchy(M, ratio)
%MP_KERNEL_CAUCHY  Family F3, interlaced Cauchy kernel A(i,j) = 1/(x_i - y_j):
%   x_i = i (i = 1..M); `ratio` points y = x_i + q/(ratio+1), q = 1..ratio, in
%   every gap  ->  N = (M-1)*ratio.  Far field: 1/(x - z_q) spans interactions with
%   sources outside the proxy circle (Cauchy integral formula).
if nargin < 2, ratio = 3; end
x = (1:M)';
q = (1:ratio)/(ratio+1);
y = reshape((x(1:end-1) + q).', [], 1);
K.m = M; K.n = numel(y);
K.x = x; K.y = y;
K.ent = @(I,J) 1 ./ (x(I) - y(J).');
K.rpos = complex(x); K.cpos = complex(y);
K.farr = @(I,zq) 1 ./ (x(I) - zq.');
K.farc = @(zq,J) 1 ./ (zq - y(J).');
K.cyclic = false; K.isreal = true;
K.name = 'cauchy_interlaced';
end
