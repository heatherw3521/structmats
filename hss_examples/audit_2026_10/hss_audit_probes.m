function T = hss_audit_probes()
%HSS_AUDIT_PROBES  Reproducers for the defects listed in HSS_Class_Audit_and_Overload_Plan.pdf.
%   T = hss_audit_probes() runs every probe and prints one line each:
%     BUG   the defect reproduces (expected on the code as of Oct 2, 2026)
%     ok    the class behaves correctly (the defect is fixed)
%     (B19 is a timing check: ok when H(1,end) grows at most ~linearly)
%   T is a struct array (id, status, detail). Each probe states the correct
%   behaviour it checks for; once a fix lands its line turns to "ok", so this
%   file doubles as a regression checklist. Runs in MATLAB and in Octave (via
%   octave_compat). Never modifies any repository file. About 1 minute.
%
%   Probe ids match the bug ids (B01, B02, ...) in the PDF.

if exist('OCTAVE_VERSION', 'builtin'), rand('seed', 1); randn('seed', 1); else, rng(1); end
ws = warning('off', 'all'); cleanup = onCleanup(@() warning(ws));
T = struct('id', {}, 'status', {}, 'detail', {});

n = 64; x = (1:n)'/n; A = 1./(x + x' + 0.5) + n*eye(n);
H = hss(A, 'blocksize', 16);

% B01 size(H) must be the matrix size
T = rep(T, 'B01', isequal(size(H), [n n]), sprintf('size(H) = [%s], expected [%d %d]', num2str(size(H)), n, n));
% B02 H(end,end) must be A(end,end)
v = tryval(@() H(end, end));
T = rep(T, 'B02', isequal(v, A(end,end)) || abs(v - A(end,end)) < 1e-12*abs(A(end,end)), sprintf('H(end,end) = %s; A(end,end) = %.6g, A(1,1) = %.6g', num2str(v), A(end,end), A(1,1)));
% B03 out-of-range / non-integer indices must error
e1 = errs(@() H(n+1, 1)); e2 = errs(@() H(0, 1)); e3 = errs(@() H(2.5, 1));
T = rep(T, 'B03', e1 && e2 && e3, sprintf('errors raised: H(n+1,1) %d, H(0,1) %d, H(2.5,1) %d', e1, e2, e3));
% B04 logical and linear indexing must follow MATLAB semantics
mask = false(n,1); mask([3 40]) = true;
s = tryval(@() H(mask, 1:3)); sz4 = size(s); ok1 = isequal(sz4, [2 3]) && norm(s - A(mask,1:3)) < 1e-10;
s = tryval(@() H(5)); ok2 = isscalar(s) && abs(s - A(5)) < 1e-10;
T = rep(T, 'B04', ok1 && ok2, sprintf('H(mask,1:3) (2 true entries) returns size %s, correct %d; H(5) correct %d', mat2str(sz4), ok1, ok2));

% B05 internal sub-blocks are not usable matrices (global indices, no root)
A5 = nested_lowrank(512, 512, 32, 5); H5 = hss(A5, 'blocksize', 32); As = A5(257:512, 257:512);
S = H5.A22;
% fixed when a piece of the tree is either refused with hss:notRoot or
% (if a subtree method is ever added) gives the right answer
[s, ~, ~, e1] = tryval(@() S(1:3, 1:3)); ok1 = strcmp(e1, 'hss:notRoot') || (isequal(size(s), [3 3]) && norm(s - As(1:3,1:3)) < 1e-10);
[y, ~, ~, e2] = tryval(@() S*ones(256,1)); ok2 = strcmp(e2, 'hss:notRoot') || (isequal(size(y), [256 1]) && norm(y - As*ones(256,1)) < 1e-10*norm(As*ones(256,1)));
T = rep(T, 'B05', ok1 && ok2, sprintf('H.A22(1:3,1:3): error "%s"; H.A22*v: error "%s" (want hss:notRoot or a correct result)', e1, e2));

% B06 threshold-mode ID silently caps ranks near 60-120 (inter_decompv3)
A6 = nested_lowrank(1024, 1024, 256, 70); H6 = hss(A6, 'blocksize', 256, 'tol', 1e-12);
w = randn(1024, 1); err = norm(H6*w - A6*w)/norm(A6*w);
T = rep(T, 'B06', err < 1e-10, sprintf('exact HSS, coupling rank 70 (leaf block-row rank 140): matvec rel. error %.2e, leaf rank found %d', err, size(H6.A11.A12.Z, 2)));

% B07 products grow ranks; the solver then rejects them (no recompression)
n7 = 256; x7 = (1:n7)'/n7; A7 = 1./(x7 + x7' + 0.5) + n7*eye(n7);
H7 = hss(A7, 'blocksize', 16); P3 = (H7*H7)*H7;
ok = ~errs(@() P3 \ randn(n7, 1));
T = rep(T, 'B07', ok, sprintf('(H*H*H)\\b: leaf rank %d > 16 leaf rows -> improperRanks error', size(P3.A11.A11.A11.A12.Z, 2)));

% B08 cutrule is not validated: k-1 split gives a wrong-size matrix
Hc = tryval(@() hss(randn(64), 'blocksize', 8, 'cutrule', @(k) k-1));
if isa(Hc, 'hss'), sz = size(tryval(@() full(Hc))); ok = isequal(sz, [64 64]); d = sprintf('full(hss(randn(64), cutrule=@(k) k-1)) is %s', mat2str(sz));
else, ok = true; d = 'rejected with an error (ok)'; end
T = rep(T, 'B08', ok, d);
% B09 sizeA option is accepted even when it contradicts size(A)
[Hs, ~, ~, eid] = tryval(@() hss(randn(64), 'blocksize', 8, 'sizeA', [32 32]));
ok = ~isa(Hs, 'hss') || isequal(size(Hs), [64 64]);
if isa(Hs, 'hss'), d = sprintf('hss(randn(64), sizeA=[32 32]) accepted; size(H) = %s', mat2str(size(Hs))); else, d = sprintf('hss(randn(64), sizeA=[32 32]) rejected: %s', eid); end
T = rep(T, 'B09', ok, d);

% B10 singular square H: silent wrong answer
A10 = nested_lowrank(256, 256, 32, 5); A10(:,7) = 0; A10(7,:) = 0;
H10 = hss(A10, 'blocksize', 32); b = randn(256, 1);
warnon(); lastwarn('');
[x10, threw, ~, eid] = tryval(@() H10 \ b); [wmsg, wid] = lastwarn(); warning('off', 'all');
named = (threw && strncmp(eid, 'hss', 3)) || strncmp(wid, 'hss', 3);
T = rep(T, 'B10', named, sprintf('A(7,:) = 0: residual %.1e; last warning [%s] "%s" (want a named hss: warning or error)', relres(A10, x10, b), wid, wmsg));
% B11 rank-deficient wide H with an inconsistent b: silent
A11 = nested_lowrank(256, 512, 32, 5); A11(200,:) = A11(10,:);
H11 = hss(A11, 'blocksize', 32); b = randn(256, 1);
warnon(); lastwarn('');
[x11, threw, ~, eid] = tryval(@() H11 \ b); wmsg = lastwarn(); warning('off', 'all');
flagged = threw || ~isempty(wmsg);
if threw, d = sprintf('two equal rows, random b: error %s', eid); else, d = sprintf('two equal rows, random b: residual %.1e, warning "%s"', relres(A11, x11, b), wmsg); end
T = rep(T, 'B11', flagged, d);
% B12 weights: singular leaf block / NaN not rejected by the check
nl = 2^H11.levelcount; Lw = repmat({eye(64)}, 1, nl); Lw{3}(5,:) = 0;
[~, ok, ~, eid] = tryval(@() minnorm(hss(nested_lowrank(256,512,32,5),'blocksize',32), randn(256,1), 'Weight', Lw));
T = rep(T, 'B12', ok, sprintf('Weight cell with a singular block: error raised %d (%s)', ok, eid));

% B13 NaN/Inf entries are accepted by the constructor
An = randn(128); An(5,100) = NaN;
[~, ok, ~, eid] = tryval(@() hss(An, 'blocksize', 16));
T = rep(T, 'B13', ok, sprintf('hss(A) with a NaN entry: error raised %d (%s)', ok, eid));
% B14 sparse input
Ssp = sprandn(256, 256, 0.02) + 10*speye(256);
[Hsp, threw] = tryval(@() hss(Ssp, 'blocksize', 32));
if threw, ok = false; d = 'hss(sparse A) errors';
else, w = randn(256,1); e = norm(Hsp*w - Ssp*w)/norm(Ssp*w); ok = e < 1e-10; d = sprintf('hss(sparse A): matvec rel. error %.1e (sparse pivoted QR is not rank revealing)', e); end
T = rep(T, 'B14', ok, d);

% B15 single-leaf root: Ir/Ic empty -> spy fails
Hl = hss(A(1:20,1:20), 'blocksize', 16);
ok = ~errs(@() spyquiet(Hl));
T = rep(T, 'B15', ok, sprintf('leaf root has Ir = [%s], Ic = [%s]; spy(H) errors', num2str(Hl.Ir), num2str(Hl.Ic)));
% B16 1-row H with scalar b
H1 = hss(randn(1, 50), 'blocksize', 16);
ok = ~errs(@() H1 \ 3);
T = rep(T, 'B16', ok, 'H\b with a 1-row H and scalar b hits the "scalar" branch and errors');
% B17 operators not overloaded
oks = [~errs(@() 2*H), ~errs(@() H*2), ~errs(@() -H), ~errs(@() H + H), ~errs(@() randn(1,n)*H)];
T = rep(T, 'B17a', all(oks), sprintf('work: 2*H %d, H*2 %d, -H %d, H+H %d, x''*H %d', oks));
ok = ~errs(@() randn(3,n)/H);
T = rep(T, 'B17b', ok, sprintf('B/H works: %d', ok));
% B18 dimension errors come from deep inside hss_matvec
[~, ~, msg] = tryval(@() H*randn(n+5, 1));
T = rep(T, 'B18', ~isempty(strfind(lower(msg), 'dimension')) || ~isempty(strfind(lower(msg), 'inner')), sprintf('H*v, wrong length: "%s"', msg));

% B19 extract builds whole off-diagonal blocks (O(n^2) per call)
tt = zeros(1, 2); nn = [2048 8192];
for k = 1:2
  xk = (1:nn(k))'/nn(k); Hk = hss(1./(xk + xk' + 0.5) + nn(k)*eye(nn(k)), 'blocksize', 64);
  tic; for r = 1:3, Hk(1, nn(k)); end; tt(k) = toc/3;
end
T = rep(T, 'B19', tt(2)/tt(1) <= 6, sprintf('H(1,end): %.4fs at n=2048, %.4fs at n=8192 (ratio %.1f; 4 = linear, 16 = quadratic)', tt(1), tt(2), tt(2)/tt(1)));

% B20 portability: sparsesign is a MEX file for Apple silicon only
d = fileparts(which('hssutil.inter_decompv3'));
hasmex = ~isempty(dir(fullfile(d, ['sparsesign.' mexext()]))) ;
hasm = exist(fullfile(d, 'sparsesign.m'), 'file') == 2;
T = rep(T, 'B20', hasm, sprintf('+hssutil has sparsesign.%s (this platform): %d; portable sparsesign.m fallback: %d', mexext(), hasmex, hasm));
% B21 two classes named hss on the path (hm-toolbox, modifed_HMtoolbox_files2019)
if ~exist('OCTAVE_VERSION', 'builtin')
  w = which('hss', '-all'); nh = sum(~cellfun(@isempty, regexp(w, '@hss[\\/]hss\.m$')));
  T = rep(T, 'B21', nh <= 1, sprintf('%d hss class folders on the path', nh));
end

fprintf('\n%-4s  %-5s  %s\n', 'id', 'state', 'detail');
for k = 1:numel(T), fprintf('%-4s  %-5s  %s\n', T(k).id, T(k).status, T(k).detail); end
fprintf('\n%d BUG, %d ok, %d info\n', sum(strcmp({T.status}, 'BUG')), sum(strcmp({T.status}, 'ok')), sum(strcmp({T.status}, 'info')));
end

% ---------------------------------------------------------------------------
function T = rep(T, id, ok, detail)
if isnan(ok), st = 'info'; elseif ok, st = 'ok'; else, st = 'BUG'; end
T(end+1) = struct('id', id, 'status', st, 'detail', detail);
end
function [v, threw, msg, id] = tryval(f)
threw = false; msg = ''; id = '';
try, v = f(); catch err, v = NaN; threw = true; msg = err.message; id = err.identifier; end
end
function warnon()
% all warnings on, except Octave's notices about the dictionary shim / parsing
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
