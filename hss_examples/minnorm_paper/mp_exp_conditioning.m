function S = mp_exp_conditioning(opts)
%MP_EXP_CONDITIONING  E2: accuracy vs conditioning.
%   Compares H\b with three other ways of computing the minimum-norm
%   solution: the dense QR solution (QR of H', the reference), the dense
%   normal equations (Cholesky of H*H') and CGNE (pcg on H*H' with HSS
%   products; flag ~= 0 means pcg did not converge).
%   Families: (a) F2 with graded columns, H*diag(10.^(-alpha*u_j)), u_j
%             uniform on [0,1] (exactly HSS for every alpha);
%             (b) F4 'valid' Gaussian convolution, stride 1, growing sigma.
%   Forward error = relative difference from the dense QR solution (both are
%   accurate to about kappa*eps, so this cannot drop below that). Backward
%   error = ||b - Hx|| / (||H|| ||x|| + ||b||), which needs no reference.
%   Writes results/conditioning.txt.
if nargin < 1, opts = struct(); end
quick = isfield(opts, 'quick') && opts.quick;
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
cgmax = 2000; if quick, cgmax = 300; end
rng(7);
S.graded = []; S.blur = [];
G0 = gen_random(128, 256, 6, 16); u = rand(256,1);
alphas = 0:2:24; sigs = 0.6:0.2:2.2;
if quick, alphas = 0:8:24; sigs = [0.6 1.2 1.8 2.2]; end
for alpha = alphas
    H = hss_from_generators(gen_scale(G0, [], 10.^(-alpha*u)));
    S.graded = [S.graded, onecond(H, alpha, cgmax)];
end
for sig = sigs
    K = kernel_conv(300, 1, 'gauss', sig, ceil(5*sig));
    [L, rb, cb] = cluster_tree(K.m, K.n, 16);
    H = hss_from_generators(hss_from_kernel(K, L, rb, cb, 1e-15));
    S.blur = [S.blur, onecond(H, sig, cgmax)];
end
T = [S.graded, S.blur];
fam = [repmat({'F2 graded'}, 1, numel(S.graded)), repmat({'F4 blur'}, 1, numel(S.blur))];
for k = 1:numel(T), T(k).family = fam{k}; end
tag = ''; if quick, tag = '_quick'; end
mp_write_table(fullfile(here, 'results', ['conditioning' tag '.txt']), T, ...
    {'family', 'param', 'kappa', 'ulv_fwd', 'ulv_bwd', 'qr_bwd', 'ne_fwd', 'ne_bwd', 'cg_fwd', 'cg_bwd', 'cg_its', 'cg_flag'}, ...
    'E2 conditioning: param = alpha (F2 graded) or sigma (F4 blur); fwd = difference from dense QR');
end

function r = onecond(H, p, cgmax)
A = full(H); b = randn(size(A,1),1); s = svd(A); nA = s(1);
xq = minnorm_dense(A, b);
bw = @(x) norm(b - A*x)/(nA*norm(x) + norm(b));
x = H \ b;
r.family = ''; r.param = p; r.kappa = s(1)/s(end);
r.ulv_fwd = relative_error(x, xq); r.ulv_bwd = bw(x); r.qr_bwd = bw(xq);
[Rc, flag] = chol(A*A');
if flag == 0
    xn = A'*(Rc \ (Rc' \ b)); r.ne_fwd = relative_error(xn, xq); r.ne_bwd = bw(xn);
else
    r.ne_fwd = NaN; r.ne_bwd = NaN;          % Cholesky of H*H' failed (kappa^2 > 1/eps)
end
Ht = H';
[y, cgflag, ~, it] = pcg(@(v) H*(Ht*v), b, 1e-14, cgmax);
xc = Ht*y; r.cg_fwd = relative_error(xc, xq); r.cg_bwd = bw(xc); r.cg_its = it; r.cg_flag = cgflag;
if flag == 0, nes = sprintf('%.1e', r.ne_fwd); else, nes = 'chol fails'; end
fprintf(['kappa=%.1e  H\\b fwd %.1e bwd %.1e | NE fwd %s | ' ...
         'CGNE fwd %.1e (flag %d, it %d) | QR bwd %.1e\n'], ...
    r.kappa, r.ulv_fwd, r.ulv_bwd, nes, r.cg_fwd, cgflag, it, r.qr_bwd);
end
