function T = hss_edge_cases()
%HSS_EDGE_CASES  Edge-case behaviour checks for @hss.
%   T = hss_edge_cases() runs every check and prints one line each:
%     ok     the class behaves as stated in the comment above the check
%     open   known limitation, not handled yet
%   T is a struct array (id, status, detail). Runs in MATLAB and in Octave
%   (via minnorm_paper/octave_compat). About 1 minute.

if exist('OCTAVE_VERSION', 'builtin'), rand('seed', 1); randn('seed', 1); else, rng(1); end
ws = warning('off', 'all'); cleanup = onCleanup(@() warning(ws));
T = struct('id', {}, 'status', {}, 'detail', {});

n = 64; x = (1:n)'/n; A = 1./(x + x' + 0.5) + n*eye(n);
H = hss(A, 'blocksize', 16);

% size(H) must be the matrix size
T = rep(T, 'size', isequal(size(H), [n n]), sprintf('size(H) = [%s], expected [%d %d]', num2str(size(H)), n, n));
% H(end,end) must be A(end,end)
v = tryval(@() H(end, end));
T = rep(T, 'end', isequal(v, A(end,end)) || abs(v - A(end,end)) < 1e-12*abs(A(end,end)), sprintf('H(end,end) = %s; A(end,end) = %.6g, A(1,1) = %.6g', num2str(v), A(end,end), A(1,1)));
% out-of-range / non-integer indices must error
e1 = errs(@() H(n+1, 1)); e2 = errs(@() H(0, 1)); e3 = errs(@() H(2.5, 1));
T = rep(T, 'range', e1 && e2 && e3, sprintf('errors raised: H(n+1,1) %d, H(0,1) %d, H(2.5,1) %d', e1, e2, e3));
% logical and linear indexing must follow MATLAB semantics
mask = false(n,1); mask([3 40]) = true;
s = tryval(@() H(mask, 1:3)); sz4 = size(s); ok1 = isequal(sz4, [2 3]) && norm(s - A(mask,1:3)) < 1e-10;
s = tryval(@() H(5)); ok2 = isscalar(s) && abs(s - A(5)) < 1e-10;
T = rep(T, 'logical', ok1 && ok2, sprintf('H(mask,1:3) (2 true entries) returns size %s, correct %d; H(5) correct %d', mat2str(sz4), ok1, ok2));

% a piece of the tree (H.A22) must be refused with hss:notRoot or give the right answer
A5 = nested_lowrank(512, 512, 32, 5); H5 = hss(A5, 'blocksize', 32); As = A5(257:512, 257:512);
S = H5.A22;
[s, ~, ~, e1] = tryval(@() S(1:3, 1:3)); ok1 = strcmp(e1, 'hss:notRoot') || (isequal(size(s), [3 3]) && norm(s - As(1:3,1:3)) < 1e-10);
[y, ~, ~, e2] = tryval(@() S*ones(256,1)); ok2 = strcmp(e2, 'hss:notRoot') || (isequal(size(y), [256 1]) && norm(y - As*ones(256,1)) < 1e-10*norm(As*ones(256,1)));
T = rep(T, 'subtree', ok1 && ok2, sprintf('H.A22(1:3,1:3): error "%s"; H.A22*v: error "%s" (want hss:notRoot or a correct result)', e1, e2));

% threshold-mode ID must find ranks above 60 (inter_decompv3)
A6 = nested_lowrank(1024, 1024, 256, 70); H6 = hss(A6, 'blocksize', 256, 'tol', 1e-12);
w = randn(1024, 1); err = norm(H6*w - A6*w)/norm(A6*w);
T = rep(T, 'idrank', err < 1e-10, sprintf('exact HSS, coupling rank 70 (leaf block-row rank 140): matvec rel. error %.2e, leaf rank found %d', err, size(H6.A11.A12.Z, 2)));

% products grow ranks; H\b on a product needs recompression first
n7 = 256; x7 = (1:n7)'/n7; A7 = 1./(x7 + x7' + 0.5) + n7*eye(n7);
H7 = hss(A7, 'blocksize', 16); P3 = (H7*H7)*H7;
ok = ~errs(@() P3 \ randn(n7, 1));
T = rep(T, 'product', ok, sprintf('(H*H*H)\\b, leaf rank %d (16 leaf rows)', size(P3.A11.A11.A11.A12.Z, 2)));

% cutrule must split properly or be rejected
Hc = tryval(@() hss(randn(64), 'blocksize', 8, 'cutrule', @(k) k-1));
if isa(Hc, 'hss'), sz = size(tryval(@() full(Hc))); ok = isequal(sz, [64 64]); d = sprintf('full(hss(randn(64), cutrule=@(k) k-1)) is %s', mat2str(sz));
else, ok = true; d = 'rejected with an error (ok)'; end
T = rep(T, 'cutrule', ok, d);
% sizeA must agree with size(A)
[Hs, ~, ~, eid] = tryval(@() hss(randn(64), 'blocksize', 8, 'sizeA', [32 32]));
ok = ~isa(Hs, 'hss') || isequal(size(Hs), [64 64]);
if isa(Hs, 'hss'), d = sprintf('hss(randn(64), sizeA=[32 32]) accepted; size(H) = %s', mat2str(size(Hs))); else, d = sprintf('hss(randn(64), sizeA=[32 32]) rejected: %s', eid); end
T = rep(T, 'sizeA', ok, d);

% singular square H must raise a named warning or error
A10 = nested_lowrank(256, 256, 32, 5); A10(:,7) = 0; A10(7,:) = 0;
H10 = hss(A10, 'blocksize', 32); b = randn(256, 1);
warnon(); lastwarn('');
[x10, threw, ~, eid] = tryval(@() H10 \ b); [wmsg, wid] = lastwarn(); warning('off', 'all');
named = (threw && strncmp(eid, 'hss', 3)) || strncmp(wid, 'hss', 3);
T = rep(T, 'singular', named, sprintf('A(7,:) = 0: residual %.1e; last warning [%s] "%s" (want a named hss: warning or error)', relres(A10, x10, b), wid, wmsg));
% rank-deficient wide H with an inconsistent b must warn or error
A11 = nested_lowrank(256, 512, 32, 5); A11(200,:) = A11(10,:);
H11 = hss(A11, 'blocksize', 32); b = randn(256, 1);
warnon(); lastwarn('');
[x11, threw, ~, eid] = tryval(@() H11 \ b); wmsg = lastwarn(); warning('off', 'all');
flagged = threw || ~isempty(wmsg);
if threw, d = sprintf('two equal rows, random b: error %s', eid); else, d = sprintf('two equal rows, random b: residual %.1e, warning "%s"', relres(A11, x11, b), wmsg); end
T = rep(T, 'inconsist', flagged, d);
% weights: a singular leaf block must be rejected
nl = 2^H11.levelcount; Lw = repmat({eye(64)}, 1, nl); Lw{3}(5,:) = 0;
[~, ok, ~, eid] = tryval(@() minnorm(hss(nested_lowrank(256,512,32,5),'blocksize',32), randn(256,1), 'Weight', Lw));
T = rep(T, 'weight', ok, sprintf('Weight cell with a singular block: error raised %d (%s)', ok, eid));

% NaN/Inf entries must be rejected by the constructor
An = randn(128); An(5,100) = NaN;
[~, ok, ~, eid] = tryval(@() hss(An, 'blocksize', 16));
T = rep(T, 'nonfinite', ok, sprintf('hss(A) with a NaN entry: error raised %d (%s)', ok, eid));
% sparse input must compress accurately
Ssp = sprandn(256, 256, 0.02) + 10*speye(256);
[Hsp, threw] = tryval(@() hss(Ssp, 'blocksize', 32));
if threw, ok = false; d = 'hss(sparse A) errors';
else, w = randn(256,1); e = norm(Hsp*w - Ssp*w)/norm(Ssp*w); ok = e < 1e-10; d = sprintf('hss(sparse A): matvec rel. error %.1e (sparse pivoted QR is not rank revealing)', e); end
T = rep(T, 'sparse', ok, d);

% spy must work on a single-leaf H
Hl = hss(A(1:20,1:20), 'blocksize', 16);
ok = ~errs(@() spyquiet(Hl));
T = rep(T, 'leafroot', ok, sprintf('spy(H) on a single-leaf H (Ir = [%s], Ic = [%s])', num2str(Hl.Ir), num2str(Hl.Ic)));
% H\b with a 1-row H and a scalar b
H1 = hss(randn(1, 50), 'blocksize', 16);
ok = ~errs(@() H1 \ 3);
T = rep(T, 'onerow', ok, 'H\b with a 1-row H and a scalar b');
% scalar products, unary minus, plus, row vector * H, B/H
oks = [~errs(@() 2*H), ~errs(@() H*2), ~errs(@() -H), ~errs(@() H + H), ~errs(@() randn(1,n)*H)];
T = rep(T, 'ops', all(oks), sprintf('2*H %d, H*2 %d, -H %d, H+H %d, x''*H %d', oks));
ok = ~errs(@() randn(3,n)/H);
T = rep(T, 'mrdivide', ok, sprintf('B/H works: %d', ok));
% a wrong-length vector in H*v must give a clear dimension error
[~, ~, msg] = tryval(@() H*randn(n+5, 1));
T = rep(T, 'dims', ~isempty(strfind(lower(msg), 'dimension')) || ~isempty(strfind(lower(msg), 'inner')), sprintf('H*v, wrong length: "%s"', msg));

% H(1,end) must cost about O(n), not O(n^2)
tt = zeros(1, 2); nn = [2048 8192];
for k = 1:2
  xk = (1:nn(k))'/nn(k); Hk = hss(1./(xk + xk' + 0.5) + nn(k)*eye(nn(k)), 'blocksize', 64);
  tic; for r = 1:3, Hk(1, nn(k)); end; tt(k) = toc/3;
end
T = rep(T, 'extract', tt(2)/tt(1) <= 6, sprintf('H(1,end): %.4fs at n=2048, %.4fs at n=8192 (ratio %.1f; 4 = linear, 16 = quadratic)', tt(1), tt(2), tt(2)/tt(1)));

% a portable sparsesign.m fallback must exist next to the MEX file
d = fileparts(which('hssutil.inter_decompv3'));
hasmex = ~isempty(dir(fullfile(d, ['sparsesign.' mexext()]))) ;
hasm = exist(fullfile(d, 'sparsesign.m'), 'file') == 2;
T = rep(T, 'sparsesign', hasm, sprintf('+hssutil has sparsesign.%s (this platform): %d; portable sparsesign.m fallback: %d', mexext(), hasmex, hasm));
% only one hss class on the path
if ~exist('OCTAVE_VERSION', 'builtin')
  w = which('hss', '-all'); nh = sum(~cellfun(@isempty, regexp(w, '@hss[\\/]hss\.m$')));
  T = rep(T, 'classes', nh <= 1, sprintf('%d hss class folders on the path', nh));
end

fprintf('\n%-10s  %-5s  %s\n', 'check', 'state', 'detail');
for k = 1:numel(T), fprintf('%-10s  %-5s  %s\n', T(k).id, T(k).status, T(k).detail); end
fprintf('\n%d ok, %d open\n', sum(strcmp({T.status}, 'ok')), sum(strcmp({T.status}, 'open')));
end

% ---------------------------------------------------------------------------
function T = rep(T, id, ok, detail)
if ok, st = 'ok'; else, st = 'open'; end
T(end+1) = struct('id', id, 'status', st, 'detail', detail);
end
function [v, threw, msg, id] = tryval(f)
threw = false; msg = ''; id = '';
try, v = f(); catch err, v = NaN; threw = true; msg = err.message; id = err.identifier; end
end
function warnon()
% all warnings on, except Octave's notices about the dictionary stand-in / parsing
warning('on', 'all');
if exist('OCTAVE_VERSION', 'builtin')
  warning('off', 'Octave:language-extension'); warning('off', 'Octave:missing-semicolon');
end
end
function t = errs(f)
t = false;
try, f(); catch, t = true; end
end
function r = relres(A, x, b)
if any(~isfinite(x(:))) || numel(x) ~= size(A,2), r = NaN; else, r = norm(A*x - b)/norm(b); end
end
function spyquiet(H)
f = figure('Visible', 'off'); c = onCleanup(@() close(f)); spy(H);
end
function A = nested_lowrank(m, n, minsize, k)
% exact HSS test matrix: rank-k coupling at every level (as in test_hss.m)
if m <= minsize || n <= minsize
  [Q1, ~] = qr(randn(m)); [Q2, ~] = qr(randn(n)); A = Q1 * eye(m, n) * Q2'; return
end
m1 = ceil(m/2); m2 = m - m1; n1 = ceil(n/2); n2 = n - n1;
A = [nested_lowrank(m1, n1, minsize, k), randn(m1, k)*randn(k, n2); ...
     randn(m2, k)*randn(k, n1), nested_lowrank(m2, n2, minsize, k)];
end
