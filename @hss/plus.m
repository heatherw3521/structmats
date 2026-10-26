function H = plus(H1, H2)
%PLUS  H1 + H2 for two HSS matrices of the same size built on the same
%   cluster tree (same blocksize and cut rule). The ranks of the sum add, so
%   the result is recompressed (see compress, tol = 1e-12). Adding a dense
%   matrix or a scalar to an HSS matrix is not supported.
if ~(isa(H1, 'hss') && isa(H2, 'hss'))
    error('hss:plus:operand', ['H + X is only supported for two HSS matrices; ' ...
        'for a dense X use full(H) + X.']);
end
hss_assertroot(H1, 'H1 + H2'); hss_assertroot(H2, 'H1 + H2');
if ~isequal(size(H1), size(H2))
    error('hss:plus:size', 'Matrix dimensions must agree (%s vs %s).', ...
        mat2str(size(H1)), mat2str(size(H2)));
end
if H1.levelcount == 0 && H2.levelcount == 0
    H = H1;
    H.D = H1.D + H2.D;
    return
end
G = hss_gen_plus(hss_to_gen(H1), hss_to_gen(H2));
G = hss_gen_compress(G, 1e-12);
G.blocksize = H1.blocksize;
H = hss_from_gen(G);
end
