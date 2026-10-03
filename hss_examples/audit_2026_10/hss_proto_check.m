function T = hss_proto_check()
%HSS_PROTO_CHECK  Accuracy checks for the prototype algorithms in ./proto that
%   back the proposed overloads in HSS_Class_Audit_and_Overload_Plan.pdf:
%     hssp_scale, hssp_plus    s*H, H1 + H2 on generators (plan: mtimes/plus)
%     hssp_orth, hssp_compress orthonormalize / recompress (plan: compress)
%     hssp_extract             H(I,J) without building blocks (plan: subsref)
%     hssp_matvec              H*X with per-level cell arrays (plan: flat storage)
%     hssp_fro                 exact Frobenius norm (plan: norm)
%     hssp_ulv_adjoint         adjoint of the stored solve (plan: tall H\b, H'\b, B/H)
%     hssp_logdet              log|det| and phase from the ULV factors (plan: det/logdet)
%   These are prototypes on the generator representation, NOT class methods;
%   none of them touches @hss. Runs in MATLAB and Octave.
here = fileparts(mfilename('fullpath')); addpath(fullfile(here, 'proto'));
if exist('OCTAVE_VERSION', 'builtin'), rand('seed', 7); randn('seed', 7); else, rng(7); end
ws = warning('off', 'all'); c = onCleanup(@() warning(ws));
T = struct('check', {}, 'value', {}, 'limit', {}, 'pass', {});

n = 512; x = (1:n)'/n; A = 1./(x + x' + 0.5) + n*eye(n);
H = hss(A, 'blocksize', 16); G = hssp_to_gen(H); FH = full(H);
B = 1./(x + x' + 1.5).^2; HB = hss(B, 'blocksize', 16); FB = full(HB);

% scale / plus
s = -3 + 1i;
F = full(hssp_from_gen(hssp_plus(hssp_scale(G, 2), hssp_scale(hssp_to_gen(HB), s))));
T = add(T, '2H + (-3+i)HB vs dense', rel(F, 2*FH + s*FB), 1e-13);
% orth is exact and gives orthonormal bases
Go = hssp_orth(G); e = 0;
for i = 1:numel(Go.U{end}), e = max(e, norm(Go.U{end}{i}'*Go.U{end}{i} - eye(size(Go.U{end}{i}, 2)))); end
T = add(T, 'orth: change of the matrix', rel(full(hssp_from_gen(Go)), FH), 1e-13);
T = add(T, 'orth: max ||U''U - I|| (leaves)', e, 1e-13);
% compress the cube: ranks drop back, solver accepts it
Q = (H*H)*H; GQ = hssp_to_gen(Q); GQc = hssp_compress(GQ, 1e-12); Qc = hssp_from_gen(GQc);
A3 = FH^3; b = randn(n, 1);
T = add(T, 'compress(H^3): rel. error', rel(full(Qc), A3), 1e-11);
T = add(T, 'compress(H^3): leaf rank (was 18)', size(GQc.U{end}{1}, 2), 8);
T = add(T, 'compress(H^3) \ b vs dense', rel(Qc \ b, A3 \ b), 1e-10);
% extract
I = [n 1 7 7 260]; J = [3 n n-1 1];
T = add(T, 'extract H(I,J) (repeats, any order)', rel(hssp_extract(G, I, J), FH(I, J)), 1e-13);
T = add(T, 'extract full row', rel(hssp_extract(G, 1, 1:n), FH(1, :)), 1e-13);
oor = false; try, hssp_extract(G, n+1, 1); catch, oor = true; end
T = add(T, 'extract: out-of-range index accepted (0 = no)', double(~oor), 0);
% matvec with cell arrays vs the class (accuracy + timing, printed)
Gm = hssp_to_gen(hss(nested_lowrank(4096, 4096, 64, 8), 'blocksize', 64)); Hm = hssp_from_gen(Gm);
xm = randn(4096, 2); tic; y1 = Hm*xm; t1 = toc; tic; y2 = hssp_matvec(Gm, xm); t2 = toc;
T = add(T, 'cell-array matvec vs H*x', rel(y2, y1), 1e-13);
fprintf('timing n = 4096, 2 columns: H*x %.4fs, hssp_matvec %.4fs\n', t1, t2);
% Frobenius norm
fr = sqrt(sum(abs(FH(:)).^2));
T = add(T, 'fro norm vs sum of squares', abs(hssp_fro(G) - fr)/fr, 1e-13);
% least squares for tall H via the adjoint of the min-norm solve of H'
At = nested_lowrank(1024, 512, 32, 6) + 1i*nested_lowrank(1024, 512, 32, 6);
Htall = hss(At, 'blocksize', 32); Hw = Htall'; Hw \ zeros(size(Hw, 1), 1);   % factor H' (fills its cache)
bt = randn(1024, 2);
xls = hssp_ulv_adjoint(Hw.factorcache.ulv, bt);
T = add(T, 'tall LS via adjoint vs pinv', rel(xls, pinv(At)*bt), 1e-11);
% H'\b reusing the factors of a square H
As = nested_lowrank(512, 512, 32, 6); Hs = hss(As, 'blocksize', 32); Hs \ zeros(512, 1);
xa = hssp_ulv_adjoint(Hs.factorcache.ulv, b);
T = add(T, 'square H''\b from H''s factors', rel(xa, As' \ b), 1e-11);
% log-determinant and phase
[ld, ph] = hssp_logdet(Hs.factorcache.ulv);
[~, Ud, Pd] = lu(As); ldr = sum(log(abs(diag(Ud)))); phr = det(Pd)*prod(sign(diag(Ud)));
T = add(T, 'logdet abs. error', abs(ld - ldr), 1e-9);
T = add(T, 'det phase error', abs(ph - phr), 1e-12);
Ac = nested_lowrank(512, 512, 32, 6) + 1i*nested_lowrank(512, 512, 32, 6); Hc = hss(Ac, 'blocksize', 32); Hc \ zeros(512, 1);
[ld, ph] = hssp_logdet(Hc.factorcache.ulv);
[~, Ud, Pd] = lu(Ac); ldr = sum(log(abs(diag(Ud)))); phr = det(Pd)*prod(diag(Ud)./abs(diag(Ud)));
T = add(T, 'logdet abs. error (complex)', abs(ld - ldr), 1e-9);
T = add(T, 'det phase error (complex)', abs(ph - phr), 1e-10);

fprintf('\n%-44s %10s %8s  %s\n', 'check', 'value', 'limit', 'result');
for k = 1:numel(T)
  r = 'PASS'; if ~T(k).pass, r = 'FAIL'; end
  fprintf('%-44s %10.2e %8.0e  %s\n', T(k).check, T(k).value, T(k).limit, r);
end
fprintf('\n%d of %d prototype checks pass\n', sum([T.pass]), numel(T));
end

function T = add(T, name, v, lim)
T(end+1) = struct('check', name, 'value', v, 'limit', lim, 'pass', v <= lim);
end
function r = rel(a, b)
r = norm(a(:) - b(:)) / norm(b(:));
end
function A = nested_lowrank(m, n, minsize, k)
if m <= minsize || n <= minsize
  [Q1, ~] = qr(randn(m)); [Q2, ~] = qr(randn(n)); A = Q1 * eye(m, n) * Q2'; return
end
m1 = ceil(m/2); m2 = m - m1; n1 = ceil(n/2); n2 = n - n1;
A = [nested_lowrank(m1, n1, minsize, k), randn(m1, k)*randn(k, n2); ...
     randn(m2, k)*randn(k, n1), nested_lowrank(m2, n2, minsize, k)];
end
