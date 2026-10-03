function K = mp_kernel_conv(n, stride, kind, sigma, w, mode)
%MP_KERNEL_CONV  Family F4, decimated 1-D convolution (deconvolution operator)
%   y_i = sum_j g(s*i + o - j) x_j,  i = 0..m-1, j = 0..n-1,
%   g(d) = exp(-d^2/(2 sigma^2))                       (kind = 'gauss')
%        = (1 - d^2/sigma^2) exp(-d^2/(2 sigma^2))     (kind = 'ricker')
%   truncated to |d| <= w (default w = ceil(4 sigma)) -> exactly banded.
%   mode 'valid': o = w, m = floor((n-1-2w)/s)+1 (stencil inside the signal)
%        'same' : o = 0, m = floor((n-1)/s)+1    (zero boundary)
%   No far field (K.farr = K.farc = []).
if nargin < 5 || isempty(w), w = ceil(4*sigma); end
if nargin < 6, mode = 'valid'; end
d = (-w:w)';
switch kind
    case 'gauss',  g = exp(-d.^2/(2*sigma^2));
    case 'ricker', g = (1 - (d/sigma).^2) .* exp(-d.^2/(2*sigma^2));
end
if strcmp(mode, 'valid'), o = w; m = floor((n-1-2*w)/stride) + 1;
else, o = 0; m = floor((n-1)/stride) + 1; end
rp = stride*(0:m-1)' + o; cp = (0:n-1)';
K.m = m; K.n = n; K.g = g; K.w = w; K.stride = stride; K.offset = o;
K.ent = @(I,J) convent(rp(I), cp(J), g, w);
K.rpos = complex(rp); K.cpos = complex(cp);
K.farr = []; K.farc = [];
K.cyclic = false; K.isreal = true;
K.name = ['conv_' kind];
end

function A = convent(r, c, g, w)
D = r - c.';
A = zeros(size(D));
mask = abs(D) <= w;
A(mask) = g(D(mask) + w + 1);
end
