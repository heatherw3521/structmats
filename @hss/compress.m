function H = compress(H, tol, mode, scale)
%COMPRESS  Recompress an HSS matrix: lower the off-diagonal ranks as far as
%   the tolerance allows. The bases are made orthonormal first, so the
%   singular values below are those of the off-diagonal block rows and
%   block columns of H itself.
%
%   H = compress(H) or compress(H, tol): relative to each block (default,
%       tol = 1e-12). In every block row and block column, singular values
%       below tol times the largest one of that block are dropped.
%   H = compress(H, tol, 'absolute'): relative to the whole matrix. Singular
%       values below tol*||H|| are dropped (||H|| estimated by normest), so
%       blocks that are small compared with H lose rank too, and the error
%       is bounded relative to ||H|| rather than block by block.
%   H = compress(H, tol, 'absolute', scale): singular values below
%       tol*scale are dropped. Use the scale of the inputs when the result
%       is much smaller than they are, e.g. after cancellation:
%       compress(H1 - H2, 1e-12, 'absolute', max(normest(H1), normest(H2))).
%
%   H1*H2 and H1 + H2 call compress(., 1e-12) automatically (their ranks
%   add); call it yourself after building H with a fixed rank or a loose
%   tolerance. Needs every leaf at the same depth (true for everything hss
%   builds).
if nargin < 2 || isempty(tol), tol = 1e-12; end
if nargin < 3 || isempty(mode), mode = 'relative'; end
hss_assertroot(H, 'compress');
if ~(isscalar(tol) && isreal(tol) && tol >= 0 && tol < 1)
    error('hss:compress:tol', 'tol must be a real scalar in [0, 1).');
end
if ~(strcmp(mode, 'relative') || strcmp(mode, 'absolute'))
    error('hss:compress:mode', 'mode must be ''relative'' or ''absolute''.');
end
if H.levelcount == 0
    return
end
atol = [];
if strcmp(mode, 'absolute')
    if nargin < 4 || isempty(scale)
        scale = normest(H);
    end
    if ~(isscalar(scale) && isreal(scale) && scale >= 0)
        error('hss:compress:scale', 'scale must be a nonnegative real scalar.');
    end
    atol = tol * scale;
end
G = hss_gen_compress(hss_to_gen(H), tol, atol);
G.blocksize = H.blocksize;
H = hss_from_gen(G);
end
