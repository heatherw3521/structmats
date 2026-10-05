classdef test_hss < matlab.unittest.TestCase
% Unit tests for the HSS matrix class.
%
% Covers:
%   1. Construction / assembly
%   2. Matrix-vector multiply    — via overloaded  H * x
%   3. Compression / rank        — tested indirectly through HSS structure
%                                  (inter_decompv3 is private; not called directly)
%   4. Solve / factorisation     — via overloaded  H \ b
%   5. Rectangular matrices      (wide solves via the minimum-norm ULV solver)
%   5b. Minimum-norm solve       — H \ b for wide H, stored factors,
%                                  minnorm (weighted) and tikhonov
%   5c. Solver on matrix families — families F1-F6 built from generators
%                                  (hss_examples/lib), vs dense
%   6.  full(H)
%   7.  Indexing H(I,J)
%   8.  Transpose H.' and H'
%   9.  HSS x HSS multiply H1*H2
%   10. Regressions: size/end, indexing, input checks, weight checks
%
% Needs hss_examples/lib (added automatically by the TestClassSetup method).
%
% Run everything:
%   results = runtests('test_hss');  table(results)
%
% Run one tag group:
%   runtests('test_hss', 'Tag', 'construction')
%   runtests('test_hss', 'Tag', 'matvec')
%   runtests('test_hss', 'Tag', 'compression')
%   runtests('test_hss', 'Tag', 'solve')
%   runtests('test_hss', 'Tag', 'rectangular')
%   runtests('test_hss', 'Tag', 'minnorm')
%   runtests('test_hss', 'Tag', 'families')
%   runtests('test_hss', 'Tag', 'full')
%   runtests('test_hss', 'Tag', 'extract')
%   runtests('test_hss', 'Tag', 'transpose')
%   runtests('test_hss', 'Tag', 'matmat')
%   runtests('test_hss', 'Tag', 'regressions')
%
% -----------------------------------------------------------------------
% leafLevel reference  (bsize = 16):
%   n = 16  ->  leafLevel = 0  (root is the only leaf)
%   n = 32  ->  leafLevel = 1  (one level of children, all leaves)
%   n = 64  ->  leafLevel = 2
%   n = 128 ->  leafLevel = 3
% Formula: leafLevel = floor( log(n/bsize) / log(2) )
% -----------------------------------------------------------------------

    properties (Constant)
        TIGHT_TOL = 1e-10   % near-exact paths (leaf dense solve, linearity)
        LOOSE_TOL = 1e-6    % randomised / multi-level paths
        BSIZE  = 16

        N_LEAF  = 16    % leafLevel = 0  (root is a single leaf)
        N_SMALL = 32    % leafLevel = 1
        N_MED   = 64    % leafLevel = 2
        N_LARGE = 128   % leafLevel = 3
    end

    methods (Static)

        function A = cauchy(m, n)
            % 1/(x_i + y_j) — naturally low-rank off-diagonals.
            if nargin < 2, n = m; end
            x = (1:m).' / m;
            y = (1:n).' / n + 0.5;
            A = 1 ./ (x + y.');
        end

        function A = spd(n)
            % Diagonally-dominant Cauchy — SPD, safe for direct solve.
            A = test_hss.cauchy(n) + n * eye(n);
        end

        function A = nestedLowRank(m, n, minsize, k)
            % Genuinely HSS-compressible rectangular matrix: every leaf is a
            % well-conditioned (full-rank) block, and every off-diagonal
            % coupling at every level is an exact rank-k correction. Unlike
            % cauchy(), whose diagonal blocks and global rank are NOT
            % guaranteed full/high, this is safe for the underdetermined
            % ULV solver, whose "immediately solved" triangular sub-blocks
            % require the diagonal/local blocks to be well-conditioned.
            if m <= minsize || n <= minsize
                [Q1,~] = qr(randn(m));
                [Q2,~] = qr(randn(n));
                A = Q1 * [eye(m), zeros(m, n-m)] * Q2';
                return
            end
            m1 = ceil(m/2); m2 = m - m1;
            n1 = ceil(n/2); n2 = n - n1;
            A11 = test_hss.nestedLowRank(m1, n1, minsize, k);
            A22 = test_hss.nestedLowRank(m2, n2, minsize, k);
            A = [A11, randn(m1,k)*randn(k,n2); randn(m2,k)*randn(k,n1), A22];
        end

        function H = build(A, varargin)
            % Build an HSS matrix using the default blocksize.
            % Extra args forwarded as named options, e.g. build(A, 'tol', 1e-10).
            H = hss(A, varargin{:}, blocksize = test_hss.BSIZE);
        end
        
        function p = leafcols(H)
            % column counts of the leaves, left to right
            if H.isleaf, p = size(H, 2); return; end
            p = [test_hss.leafcols(H.A11), test_hss.leafcols(H.A22)];
        end

        function q = leafrows(H)
            % row counts of the leaves, top to bottom
            if H.isleaf, q = size(H, 1); return; end
            q = [test_hss.leafrows(H.A11), test_hss.leafrows(H.A22)];
        end

        function d = factorDepth(F)
            % number of levels in a stored ULV factorization
            d = 0;
            while ~F.isroot
                d = d + 1;
                F = F.next;
            end
        end

        function Afull = dense(H)
            % Recover the dense matrix via the overloaded subsref / extract.
            Afull = H(1:size(H, 1), 1:size(H, 2));
        end
    end


    % ==================================================================
    % 1. CONSTRUCTION / ASSEMBLY
    % ==================================================================
    methods (TestClassSetup)
        function addHelperLibrary(tc)
            % hss_examples/lib: generators, kernels, dense references
            lib = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'lib');
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(lib));
        end
    end

    methods (Test, TestTags = {'construction'})

        function test_empty_constructor(tc)
            H = hss();
            tc.verifyEqual(size(H), [0 0], 'empty hss: size must be [0 0]');
            tc.verifyTrue(isempty(H), 'empty hss: isempty must be true');
            tc.verifyEmpty(H.D,     'empty hss: D must be []');
            tc.verifyEmpty(H.A11,   'empty hss: A11 must be []');
        end

        function test_size_matches_input(tc)
            n = tc.N_MED;
            H = tc.build(tc.cauchy(n));
            tc.verifyEqual(size(H), [n, n], ...
                'size(H) must equal input matrix dimensions');
        end

        function test_root_flags(tc)
            H = tc.build(tc.cauchy(tc.N_MED));
            tc.verifyTrue(logical(H.isroot),  'root must have isroot=true');
            tc.verifyFalse(logical(H.isleaf), 'root with n>blocksize must have isleaf=false');
            tc.verifyTrue(logical(H.isdiag),  'root must have isdiag=true');
        end

        function test_leaf_when_n_equals_blocksize(tc)
            % n == blocksize  ->  leafLevel = 0  ->  root is a leaf
            n = tc.N_LEAF;
            H = hss(tc.cauchy(n), blocksize = n);
            tc.verifyTrue(logical(H.isleaf), ...
                'root must be a leaf when n == blocksize (leafLevel=0)');
        end

        function test_leaf_stores_dense_block(tc)
            n = tc.N_LEAF;
            A = tc.cauchy(n);
            H = hss(A, blocksize = n);
            tc.verifyNotEmpty(H.D, 'leaf node must have a non-empty D block');
            tc.verifySize(H.D, [n, n], 'leaf D must be n x n');
        end

        function test_leaf_D_matches_input(tc)
            % Single-leaf tree: D must equal A exactly.
            n = tc.N_LEAF;
            A = tc.cauchy(n);
            H = hss(A, blocksize = n);
            tc.verifyEqual(H.D, A, 'AbsTol', 1e-15, ...
                'leaf D must equal the original matrix');
        end

        function test_children_exist_for_nonleaf(tc)
            H = tc.build(tc.cauchy(tc.N_MED));
            tc.verifyNotEmpty(H.A11, 'non-leaf root must have A11');
            tc.verifyNotEmpty(H.A22, 'non-leaf root must have A22');
            tc.verifyNotEmpty(H.A12, 'non-leaf root must have A12');
            tc.verifyNotEmpty(H.A21, 'non-leaf root must have A21');
        end

        function test_child_sizes_sum_to_parent(tc)
            n = tc.N_MED;
            H = tc.build(tc.cauchy(n));
            tc.verifyEqual(size(H.A11, 1) + size(H.A22, 1), n, ...
                'child row-sizes must sum to parent row-size');
            tc.verifyEqual(size(H.A11, 2) + size(H.A22, 2), n, ...
                'child col-sizes must sum to parent col-size');
        end

        function test_offdiag_sizes_consistent(tc)
            n = tc.N_MED;
            H = tc.build(tc.cauchy(n));
            tc.verifyEqual(size(H.A12, 1), size(H.A11, 1), ...
                'A12 row-size must match A11 row-size');
            tc.verifyEqual(size(H.A12, 2), size(H.A22, 2), ...
                'A12 col-size must match A22 col-size');
            tc.verifyEqual(size(H.A21, 1), size(H.A22, 1), ...
                'A21 row-size must match A22 row-size');
            tc.verifyEqual(size(H.A21, 2), size(H.A11, 2), ...
                'A21 col-size must match A11 col-size');
        end

        function test_index_ranges_span_full_matrix(tc)
            n = tc.N_MED;
            H = tc.build(tc.cauchy(n));
            tc.verifyEqual(H.Ir, [1, n], 'root Ir must be [1, n]');
            tc.verifyEqual(H.Ic, [1, n], 'root Ic must be [1, n]');
        end

        function test_child_index_ranges_contiguous(tc)
            n = tc.N_MED;
            H = tc.build(tc.cauchy(n));
            tc.verifyEqual(H.A11.Ir(2) + 1, H.A22.Ir(1), ...
                'A11 and A22 row ranges must be contiguous');
            tc.verifyEqual(H.A11.Ic(2) + 1, H.A22.Ic(1), ...
                'A11 and A22 col ranges must be contiguous');
        end

        function test_offdiag_nodes_not_isdiag(tc)
            H = tc.build(tc.cauchy(tc.N_MED));
            tc.verifyFalse(logical(H.A12.isdiag), 'A12 must have isdiag=false');
            tc.verifyFalse(logical(H.A21.isdiag), 'A21 must have isdiag=false');
        end

        function test_rowtreeindex_pattern(tc)
            % Left child gets odd index, right child gets even.
            H = tc.build(tc.cauchy(tc.N_MED));
            tc.verifyEqual(H.A11.rowtreeindex, 1, 'A11 rowtreeindex must be 1');
            tc.verifyEqual(H.A22.rowtreeindex, 2, 'A22 rowtreeindex must be 2');
        end

        function test_levelcount_at_root(tc)
            n = tc.N_MED;
            H = tc.build(tc.cauchy(n));
            expected = floor(log(n / tc.BSIZE) / log(2));
            tc.verifyEqual(H.levelcount, expected, ...
                'root levelcount must equal computed leafLevel');
        end

        function test_reconstruction_accuracy(tc)
            n = tc.N_MED;
            A = tc.cauchy(n);
            H = tc.build(A);
            relerr = norm(A - tc.dense(H), 'fro') / norm(A, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'dense reconstruction must match original matrix');
        end

        function test_identity_reconstructs_exactly(tc)
            n = tc.N_SMALL;
            A = eye(n);
            H = tc.build(A);
            relerr = norm(A - tc.dense(H), 'fro') / (norm(A, 'fro') + eps);
            tc.verifyLessThan(relerr, tc.TIGHT_TOL, ...
                'identity matrix must reconstruct exactly');
        end

        function test_k_option_prescribes_rank(tc)
            n = tc.N_MED;
            k_target = 4;
            H = hss(tc.cauchy(n), blocksize = tc.BSIZE, k = k_target);
            tc.verifyEqual(size(H.A11.A12.Z, 2), k_target, ...
                'k= option must fix the leaf off-diagonal rank to k');
        end

        function test_tol_option_accepted(tc)
            n = tc.N_MED;
            A = tc.cauchy(n);
            H = hss(A, blocksize = tc.BSIZE, tol = 1e-12);
            relerr = norm(A - tc.dense(H), 'fro') / norm(A, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'tol=1e-12 must still produce an accurate HSS');
        end

        function test_decomp_oldid_option_works(tc)
            % Regression: hss_constructor used `options.decomp == 'ID'`
            % (char-array ==) instead of strcmp, which errors outright
            % ("Arrays have incompatible sizes") whenever the option
            % string has a different length than 'ID' -- so 'oldID' (and
            % 'SVD') were completely unreachable regardless of which
            % branch the caller actually wanted. Fixed to strcmp.
            n = tc.N_MED;
            A = tc.cauchy(n);
            H = hss(A, blocksize = tc.BSIZE, decomp = 'oldID');
            relerr = norm(A - tc.dense(H), 'fro') / norm(A, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'decomp=oldID must still produce an accurate HSS');
        end

        function test_invalid_blocksize_errors(tc)
            % Regression: blocksize <= 0 made hss_constructor's leafLevel
            % search loop infinite -- it only breaks once floor(.../2)
            % drops below blocksize, which never happens for
            % blocksize <= 0. Now validated in the arguments block
            % instead of hanging.
            tc.verifyError(@() hss(tc.cauchy(tc.N_MED), blocksize = 0), ...
                'MATLAB:validators:mustBePositive', 'blocksize = 0 must error, not hang');
            tc.verifyError(@() hss(tc.cauchy(tc.N_MED), blocksize = -4), ...
                'MATLAB:validators:mustBePositive', 'negative blocksize must error, not hang');
        end

        function test_invalid_decomp_option_errors(tc)
            % An unrecognized decomp value must error rather than silently
            % fall through to a default.
            tc.verifyError(@() hss(tc.cauchy(tc.N_MED), blocksize = tc.BSIZE, decomp = 'bogus'), 'hss:decomp', ...
                'an unrecognized decomp option must error, not silently pick a default');
        end

        function test_svd_decomp_not_implemented(tc)
            % 'SVD' is a recognized option that is not implemented yet.
            tc.verifyError(@() hss(tc.cauchy(tc.N_MED), blocksize = tc.BSIZE, decomp = 'SVD'), ...
                'hss:decomp:notImplemented', 'decomp=SVD must say it is not implemented');
        end

        function test_custom_cutrule_controls_split_point(tc)
            % cutrule is otherwise only ever exercised incidentally (as a
            % way to force two HSS trees to disagree, in the matmat
            % incompatible-structure test) -- this checks hss() itself
            % actually honors it for where a level's split falls.
            n = tc.N_MED;
            A = tc.cauchy(n);
            H = hss(A, blocksize = tc.BSIZE, cutrule = @(k) floor(k/4));
            tc.verifyEqual(size(H.A11, 1), floor(n/4), ...
                'custom cutrule must control the row split point');
            tc.verifyEqual(size(H.A11, 2), floor(n/4), ...
                'custom cutrule must control the col split point');
        end

        function test_spy_runs_without_error(tc)
            % Regression: size(H, 22) (should be size(H, 2)) in spy.m's
            % draw_spy, in a branch only reached at H.level==cutoff+1 --
            % would throw "index exceeds array bounds" the moment that
            % branch executed. N_LARGE/BSIZE gives levelcount=3 so that
            % branch is actually reached here.
            H = tc.build(tc.cauchy(tc.N_LARGE));
            fig = figure('Visible', 'off');
            cleanupObj = onCleanup(@() close(fig)); %#ok<NASGU>
            spy(H);
        end
    end


    % ==================================================================
    % 2. MATRIX-VECTOR MULTIPLY  (overloaded  *)
    % ==================================================================
    methods (Test, TestTags = {'matvec'})

        function test_matvec_matches_dense(tc)
            n = tc.N_MED;
            A = tc.cauchy(n);
            H = tc.build(A);
            x = randn(n, 1);
            relerr = norm(H*x - A*x) / norm(A*x);
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, 'H*x must match A*x');
        end

        function test_matvec_zero_vector(tc)
            n = tc.N_MED;
            H = tc.build(tc.cauchy(n));
            tc.verifyEqual(H * zeros(n,1), zeros(n,1), 'AbsTol', 1e-14, ...
                'H * 0 must be 0');
        end

        function test_matvec_linearity(tc)
            n = tc.N_SMALL;
            A = tc.cauchy(n);
            H = tc.build(A);
            x1 = randn(n,1);  x2 = randn(n,1);  alpha = pi;
            lhs = H * (alpha*x1 + x2);
            rhs = alpha*(H*x1) + H*x2;
            tc.verifyLessThan(norm(lhs-rhs)/(norm(rhs)+eps), 1e-12, ...
                'H*(alpha*x1+x2) must equal alpha*H*x1 + H*x2');
        end

        function test_matvec_identity(tc)
            n = tc.N_SMALL;
            H = tc.build(eye(n));
            x = randn(n,1);
            tc.verifyLessThan(norm(H*x - x)/norm(x), tc.TIGHT_TOL, ...
                'I*x must equal x');
        end

        function test_matvec_multiple_columns(tc)
            n = tc.N_SMALL;
            A = tc.cauchy(n);
            H = tc.build(A);
            X = randn(n, 5);
            relerr = norm(H*X - A*X, 'fro') / norm(A*X, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'multi-column H*X must match A*X');
        end

        function test_matvec_consistent_across_blocksizes(tc)
            n = tc.N_MED;
            A = tc.cauchy(n);
            x = randn(n,1);
            y_ref = A * x;
            for bsize = [8, 16, 32]
                H = hss(A, blocksize = bsize);
                relerr = norm(H*x - y_ref) / norm(y_ref);
                tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                    sprintf('mat-vec must be accurate for blocksize=%d', bsize));
            end
        end

        function test_scalar_multiply_throws(tc)
            H = tc.build(tc.cauchy(tc.N_SMALL));
            tc.verifyError(@() H * 3.0, 'hss:mtimes:scalar', ...
                'H * scalar must throw (not yet supported)');
        end

        function test_left_multiply_throws(tc)
            H = tc.build(tc.cauchy(tc.N_SMALL));
            tc.verifyError(@() randn(1, tc.N_SMALL) * H, 'hss:mtimes:left', ...
                'x * H with a dense x must throw (not supported)');
        end
    end


    % ==================================================================
    % 3. COMPRESSION / RANK
    %
    % inter_decompv3 is private so it is not called directly.
    % These tests verify its effects through the public HSS interface:
    %   (a) off-diagonal factor dimensions and index bookkeeping
    %   (b) reconstruction accuracy at varying tolerances / prescribed k
    %   (c) mat-vec error as a proxy for compression quality
    % ==================================================================
    methods (Test, TestTags = {'compression'})

        function test_offdiag_Z_and_Y_nonempty(tc)
            H = tc.build(tc.cauchy(tc.N_MED));
            tc.verifyNotEmpty(H.A12.Z, 'A12.Z must be non-empty');
            tc.verifyNotEmpty(H.A12.Y, 'A12.Y must be non-empty');
            tc.verifyNotEmpty(H.A21.Z, 'A21.Z must be non-empty');
            tc.verifyNotEmpty(H.A21.Y, 'A21.Y must be non-empty');
        end

        function test_offdiag_rank_below_blocksize(tc)
            % Cauchy is HSS-compressible; leaf off-diagonal rank << blocksize.
            H = tc.build(tc.cauchy(tc.N_MED));
            k12 = size(H.A11.A12.Z, 2);
            tc.verifyLessThanOrEqual(k12, tc.BSIZE, ...
                'leaf off-diagonal rank must be <= blocksize for Cauchy input');
        end

        function test_offdiag_rank_much_less_than_n(tc)
            n = tc.N_MED;
            H = tc.build(tc.cauchy(n));
            tc.verifyLessThan(size(H.A12.Z, 2), n/2, ...
                'root off-diagonal rank must be << n/2');
        end

        function test_prescribed_k_propagates_to_leaves(tc)
            k = 3;
            H = hss(tc.cauchy(tc.N_MED), blocksize = tc.BSIZE, k = k);
            tc.verifyEqual(size(H.A11.A12.Z, 2), k, ...
                'prescribed k must appear in leaf off-diagonal blocks');
        end

        function test_lrcomponent_size_consistent_with_ZY(tc)
            % lrcomponent = A(lowrankrows, lowrankcols)
            % must be  size(Z,2) x size(Y,1)  at each leaf off-diagonal block.
            H = tc.build(tc.cauchy(tc.N_MED));
            Hleaf = H.A11.A12;
            tc.verifyEqual(size(Hleaf.lrcomponent, 1), size(Hleaf.Z, 2), ...
                'lrcomponent rows must equal Z cols (rank k)');
            tc.verifyEqual(size(Hleaf.lrcomponent, 2), size(Hleaf.Y, 1), ...
                'lrcomponent cols must equal Y rows (rank k)');
        end

        function test_lowrankrows_within_matrix_bounds(tc)
            n = tc.N_MED;
            H = tc.build(tc.cauchy(n));
            rows = H.A11.A12.lowrankrows;
            tc.verifyGreaterThanOrEqual(min(rows), 1, ...
                'lowrankrows must be >= 1');
            tc.verifyLessThanOrEqual(max(rows), n, ...
                'lowrankrows must be <= n');
        end

        function test_lowrankcols_within_matrix_bounds(tc)
            n = tc.N_MED;
            H = tc.build(tc.cauchy(n));
            cols = H.A11.A12.lowrankcols;
            tc.verifyGreaterThanOrEqual(min(cols), 1, ...
                'lowrankcols must be >= 1');
            tc.verifyLessThanOrEqual(max(cols), n, ...
                'lowrankcols must be <= n');
        end

        % function test_tighter_tol_rank_not_larger(tc)
        %     n = tc.N_MED;
        %     A = tc.cauchy(n);
        %     H_loose = hss(A, blocksize = tc.BSIZE, tol = 1e-4);
        %     H_tight = hss(A, blocksize = tc.BSIZE, tol = 1e-12);
        %     k_loose = size(H_loose.A11.A12.Z, 2);
        %     k_tight = size(H_tight.A11.A12.Z, 2);
        %     tc.verifyGreaterThanOrEqual(k_tight, k_loose, ...
        %         'tighter tolerance must not inflate off-diagonal rank');
        % end

        function test_reconstruction_improves_with_tighter_tol(tc)
            % Reconstruction error must be no worse with a tighter tolerance.
            n = tc.N_MED;
            A = tc.cauchy(n);
            H_loose = hss(A, blocksize = tc.BSIZE, tol = 1e-4);
            H_tight = hss(A, blocksize = tc.BSIZE, tol = 1e-12);
            err_loose = norm(A - tc.dense(H_loose), 'fro') / norm(A, 'fro');
            err_tight = norm(A - tc.dense(H_tight), 'fro') / norm(A, 'fro');
            tc.verifyLessThanOrEqual(err_tight, err_loose, ...
                'tighter tol should give equal-or-better reconstruction');
        end

        function test_matvec_error_reflects_compression_rank(tc)
            % Low prescribed rank -> larger mat-vec error vs high rank.
            n = tc.N_MED;
            A = tc.cauchy(n);
            x = randn(n, 1);
            y_ref = A * x;
            H_rich = hss(A, blocksize = tc.BSIZE, tol = 1e-12);
            H_poor = hss(A, blocksize = tc.BSIZE, k = 1);
            err_rich = norm(H_rich*x - y_ref) / norm(y_ref);
            err_poor = norm(H_poor*x - y_ref) / norm(y_ref);
            tc.verifyLessThan(err_rich, err_poor, ...
                'higher-rank HSS must have smaller mat-vec error than rank-1');
        end
    end


    % ==================================================================
    % 4. SOLVE  (overloaded  \)
    % ==================================================================
    methods (Test, TestTags = {'solve'})

        function test_solve_recovers_true_x(tc)
            n = tc.N_SMALL;
            A = tc.spd(n);
            H = tc.build(A);
            x_true = randn(n, 1);
            x = H \ (A * x_true);
            relerr = norm(x - x_true) / norm(x_true);
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'H\\b must recover x_true');
        end

        function test_solve_small_residual(tc)
            n = tc.N_SMALL;
            A = tc.spd(n);
            H = tc.build(A);
            b = randn(n, 1);
            x = H \ b;
            relerr = norm(A*x - b) / norm(b);
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                '||Ax - b|| / ||b|| must be small');
        end

        function test_solve_leaf_is_exact(tc)
            % leafLevel = 0: single leaf, falls through to D\b — near-exact.
            n = tc.N_LEAF;
            A = tc.spd(n);
            H = hss(A, blocksize = n);
            x_true = randn(n, 1);
            x = H \ (A * x_true);
            relerr = norm(x - x_true) / norm(x_true);
            tc.verifyLessThan(relerr, tc.TIGHT_TOL, ...
                'single-leaf solve must be near machine-exact');
        end

        function test_solve_identity_system(tc)
            n = tc.N_SMALL;
            H = tc.build(eye(n));
            b = randn(n, 1);
            x = H \ b;
            tc.verifyLessThan(norm(x - b) / norm(b), tc.TIGHT_TOL, ...
                'I\\b must return b');
        end

        function test_solve_zero_rhs(tc)
            n = tc.N_SMALL;
            H = tc.build(tc.spd(n));
            x = H \ zeros(n, 1);
            tc.verifyEqual(x, zeros(n, 1), 'AbsTol', 1e-14, ...
                'H\\0 must return 0');
        end

        function test_solve_does_not_mutate_rhs(tc)
            n = tc.N_SMALL;
            H = tc.build(tc.spd(n));
            b = randn(n, 1);
            b_copy = b;
            H \ b;
            tc.verifyEqual(b, b_copy, ...
                'H\\b must not modify the input b');
        end

        function test_solve_consistent_with_matvec(tc)
            % y = H*x_true  ->  H\y should give back x_true.
            n = tc.N_SMALL;
            A = tc.spd(n);
            H = tc.build(A);
            x_true = randn(n, 1);
            x = H \ (H * x_true);
            relerr = norm(x - x_true) / norm(x_true);
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'H \\ (H*x) must recover x');
        end

        function test_solve_multiple_rhs(tc)
            n = tc.N_SMALL;
            A = tc.spd(n);
            H = tc.build(A);
            X_true = randn(n, 3);
            B = A * X_true;
            for j = 1:3
                x = H \ B(:, j);
                relerr = norm(x - X_true(:,j)) / norm(X_true(:,j));
                tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                    sprintf('solve must be accurate for RHS column %d', j));
            end
        end

        function test_solve_robust_across_blocksizes(tc)
            n = tc.N_MED;
            A = tc.spd(n);
            x_true = randn(n, 1);
            b = A * x_true;
            for bsize = [8, 16, 32]
                H = hss(A, blocksize = bsize);
                x = H \ b;
                relerr = norm(x - x_true) / norm(x_true);
                tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                    sprintf('solve must be accurate for blocksize=%d', bsize));
            end
        end

        function test_mldivide_unsupported_combinations_error(tc)
            % mldivide.m has three branches beyond "hss \ dense vector"
            % that are all explicitly unimplemented (hss\hss, hss\scalar,
            % dense\hss / hss inverse) -- confirm each fails loudly rather
            % than silently misbehaving. Grouped into one test since all
            % three assert the same thing (throws) about the same file.
            n = tc.N_SMALL;
            H1 = tc.build(tc.spd(n));
            H2 = tc.build(tc.spd(n));
            tc.verifyError(@() H1 \ H2, 'hss:mldivide:hss', 'hss \\ hss must throw (not yet supported)');
            tc.verifyError(@() H1 \ 3.0, 'hss:mldivide:scalar', 'hss \\ scalar must throw (not yet supported)');
            tc.verifyError(@() tc.spd(n) \ H1, 'hss:mldivide:dense', 'dense \\ hss must throw (inverse not yet supported)');
        end
    end


    % ==================================================================
    % 5. RECTANGULAR MATRICES
    %
    % Construction and mat-vec are fully tested. Wide solves (m < n) go
    % through the minimum-norm ULV solver (@hss/private/
    % hss_ulvminnormsolve.m); its own tests are in section 5b. Tall
    % solves are not implemented yet.
    % ==================================================================
    methods (Test, TestTags = {'rectangular'})

        % ---- Construction -------------------------------------------

        function test_rect_tall_size(tc)
            m = 64; n = 32;
            H = hss(tc.cauchy(m,n), blocksize = tc.BSIZE);
            tc.verifyEqual(size(H), [m, n], 'tall HSS must record [m, n]');
        end

        function test_rect_wide_size(tc)
            m = 32; n = 64;
            H = hss(tc.cauchy(m,n), blocksize = tc.BSIZE);
            tc.verifyEqual(size(H), [m, n], 'wide HSS must record [m, n]');
        end

        function test_rect_root_flags(tc)
            H = hss(tc.cauchy(64,32), blocksize = tc.BSIZE);
            tc.verifyTrue(logical(H.isroot),  'rect root must have isroot=true');
            tc.verifyFalse(logical(H.isleaf), 'rect root must have isleaf=false');
            tc.verifyTrue(logical(H.isdiag),  'rect root must have isdiag=true');
        end

        function test_rect_child_sizes_sum_to_parent(tc)
            m = 64; n = 32;
            H = hss(tc.cauchy(m,n), blocksize = tc.BSIZE);
            tc.verifyEqual(size(H.A11, 1) + size(H.A22, 1), m, ...
                'rect child row-sizes must sum to m');
            tc.verifyEqual(size(H.A11, 2) + size(H.A22, 2), n, ...
                'rect child col-sizes must sum to n');
        end

        function test_rect_tall_reconstruction(tc)
            m = 64; n = 32;
            A = tc.cauchy(m, n);
            H = hss(A, blocksize = tc.BSIZE);
            relerr = norm(A - tc.dense(H), 'fro') / norm(A, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'tall HSS reconstruction must match original matrix');
        end

        function test_rect_wide_reconstruction(tc)
            m = 32; n = 64;
            A = tc.cauchy(m, n);
            H = hss(A, blocksize = tc.BSIZE);
            relerr = norm(A - tc.dense(H), 'fro') / norm(A, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'wide HSS reconstruction must match original matrix');
        end

        % ---- Mat-vec ------------------------------------------------

        function test_rect_tall_matvec(tc)
            m = 64; n = 32;
            A = tc.cauchy(m, n);
            H = hss(A, blocksize = tc.BSIZE);
            x = randn(n, 1);
            relerr = norm(H*x - A*x) / norm(A*x);
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, 'tall H*x must match A*x');
        end

        function test_rect_wide_matvec(tc)
            m = 32; n = 64;
            A = tc.cauchy(m, n);
            H = hss(A, blocksize = tc.BSIZE);
            x = randn(n, 1);
            relerr = norm(H*x - A*x) / norm(A*x);
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, 'wide H*x must match A*x');
        end

        function test_rect_matvec_output_size(tc)
            m = 64; n = 32;
            H = hss(tc.cauchy(m,n), blocksize = tc.BSIZE);
            y = H * randn(n, 1);
            tc.verifySize(y, [m, 1], 'output of tall H*x must be m x 1');
        end

        function test_rect_matvec_zero_vector(tc)
            m = 64; n = 32;
            H = hss(tc.cauchy(m,n), blocksize = tc.BSIZE);
            tc.verifyEqual(H * zeros(n,1), zeros(m,1), 'AbsTol', 1e-14, ...
                'tall H * 0 must be 0');
        end

        function test_rect_matvec_linearity(tc)
            m = 64; n = 32;
            A = tc.cauchy(m, n);
            H = hss(A, blocksize = tc.BSIZE);
            x1 = randn(n,1);  x2 = randn(n,1);  alpha = exp(1);
            lhs = H * (alpha*x1 + x2);
            rhs = alpha*(H*x1) + H*x2;
            tc.verifyLessThan(norm(lhs-rhs)/(norm(rhs)+eps), 1e-12, ...
                'rect H*(alpha*x1+x2) must equal alpha*H*x1 + H*x2');
        end

        % ---- Solve ----------------------------------------------------

        function test_rect_tall_solve_errors(tc)
            % Overdetermined solve (m > n) is explicitly unimplemented in
            % mldivide.m (unconditional error in that size branch) --
            % confirm it fails loudly instead of returning garbage.
            % When implemented: replace with a residual-based check like
            % test_rect_wide_solve's.
            m = 64; n = 32;
            A = tc.cauchy(m, n);
            H = hss(A, blocksize = tc.BSIZE);
            b = randn(m, 1);
            tc.verifyError(@() H \ b, 'hss:mldivide:tall', ...
                'overdetermined H\\b must throw (not yet implemented)');
        end

        function test_rect_wide_solve(tc)
            % Underdetermined solve (m < n), via hss_ulvminnormsolve.
            % Uses nestedLowRank rather than cauchy(): the solver needs H to
            % have full row rank, and cauchy() blocks can be numerically
            % rank deficient.
            m = 32; n = 64;
            rng(1);
            A = tc.nestedLowRank(m, n, tc.BSIZE, 2);
            H = hss(A, blocksize = tc.BSIZE, k = 2);
            b = randn(m, 1);
            x = H \ b;
            relerr = norm(H*x - b) / norm(b);
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'underdetermined H\\b must satisfy H*x = b');
        end

        function test_wide_solve_when_ulv_slack_short(tc)
            % A leaf pair whose slack (n_tau - m_tau = 2) is smaller than
            % its off-diagonal rank (6). The size reduction keeps
            % min(n_tau, l_tau + m_tau) columns, so no slack is needed;
            % section 5b checks that no level is merged and that x is the
            % minimum-norm solution.
            %
            % Uses tol= rather than k=: forcing rank 6 by truncation on
            % these blocks would make H itself rank deficient, so no x
            % would satisfy H*x = b.
            m = 16; n = 20; bs = 8; k = 6;
            rng(1);
            A = tc.nestedLowRank(m, n, bs, k);   % offdiag rank (6) exceeds leaf slack (n_tau-m_tau=2)
            H = hss(A, blocksize = bs, tol = 1e-10);
            b = randn(m, 1);
            x = H \ b;
            relerr = norm(H*x - b) / norm(b);
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'solve must satisfy Hx=b when a leaf''s slack is shorter than its off-diagonal rank');
        end
    end


    % ==================================================================
    % 5b. MINIMUM-NORM SOLVE  (H \ b for wide H, stored factors,
    %     minnorm and tikhonov)
    %
    % @hss/private/hss_ulvminnormsolve.m returns x = argmin ||x|| subject
    % to H*x = b (for square H, the solution). These tests check the
    % minimum-norm property itself (not just the residual), several
    % right-hand sides, that H keeps its factors between solves (and drops
    % them when it is modified), square systems, complex data, backward
    % stability on an ill-conditioned matrix, that no tree level is merged
    % when a leaf's slack is short, the named error for a rank-deficient H,
    % and weighted minimum norm / Tikhonov with leaf-block weights.
    % ==================================================================
    methods (Test, TestTags = {'minnorm', 'rectangular'})

        function test_minnorm_matches_pinv(tc)
            % tol= rather than k=: with two tree levels a leaf's row basis
            % has to carry both levels' rank-2 couplings
            m = 64; n = 128;
            rng(2);
            A = tc.nestedLowRank(m, n, tc.BSIZE, 2);
            H = hss(A, blocksize = tc.BSIZE, tol = 1e-12);
            Hf = full(H);
            b = randn(m, 1);
            x = H \ b;
            xref = pinv(Hf) * b;
            tc.verifyLessThan(norm(x - xref) / norm(xref), tc.TIGHT_TOL, ...
                'H\\b must be the minimum-norm solution pinv(H)*b');
            [Q, ~] = qr(Hf', 0);                 % orthonormal basis of range(H')
            tc.verifyLessThan(norm(x - Q*(Q'*x)) / norm(x), tc.TIGHT_TOL, ...
                'the minimum-norm solution must have no component in null(H)');
        end

        function test_minnorm_slack_short_merges_no_level(tc)
            % same matrix as test_wide_solve_when_ulv_slack_short
            m = 16; n = 20; bs = 8; k = 6;
            rng(1);
            A = tc.nestedLowRank(m, n, bs, k);   % leaf slack 2 < off-diagonal rank 6
            H = hss(A, blocksize = bs, tol = 1e-10);
            b = randn(m, 1);
            x = H \ b;
            F = H.factorcache.ulv;               % factors stored by H\b
            depth = 0;
            while ~F.isroot
                depth = depth + 1;
                F = F.next;
            end
            tc.verifyEqual(depth, H.levelcount, ...
                'the factorization must keep every tree level (no levelup merging)');
            Hf = full(H);
            xref = pinv(Hf) * b;
            tc.verifyLessThan(norm(x - xref) / norm(xref), tc.TIGHT_TOL, ...
                'slack-short leaves must still give the minimum-norm solution');
        end

        function test_minnorm_multiple_rhs(tc)
            m = 64; n = 128;
            rng(3);
            H = hss(tc.nestedLowRank(m, n, tc.BSIZE, 2), blocksize = tc.BSIZE, tol = 1e-12);
            B = randn(m, 3);
            X = H \ B;
            tc.verifyEqual(size(X), [n, 3], 'H\\B must return one column per right-hand side');
            for j = 1:3
                xj = H \ B(:, j);
                tc.verifyLessThan(norm(X(:, j) - xj) / norm(xj), tc.TIGHT_TOL, ...
                    'H\\B must solve each column as H\\b does');
            end
            Xr = mp_minnorm_dense(full(H), B);
            tc.verifyLessThan(norm(X - Xr) / norm(Xr), tc.TIGHT_TOL, ...
                'H\\B must be the minimum-norm solution of every column');
        end

        function test_backslash_reuses_stored_factors(tc)
            m = 64; n = 128;
            rng(4);
            H = hss(tc.nestedLowRank(m, n, tc.BSIZE, 2), blocksize = tc.BSIZE, tol = 1e-12);
            b = randn(m, 1);
            tc.verifyEmpty(H.factorcache.ulv, 'a new H must not hold factors yet');
            x = H \ b;
            tc.verifyNotEmpty(H.factorcache.ulv, 'H\\b must keep its factors with H');
            tc.verifyEqual(H \ b, x, 'a second H\\b must reuse the stored factors');
            tc.verifyEqual(tc.factorDepth(H.factorcache.ulv), H.levelcount, ...
                'the stored factors must hold one entry per tree level');
            % the projection onto {x : H*x = b} used by Douglas-Rachford/ADMM
            v = randn(n, 1);
            Pv = v + H \ (b - H*v);
            tc.verifyLessThan(norm(H*Pv - b) / norm(b), tc.TIGHT_TOL, ...
                'the projection P(v) = v + H^+(b - H*v) must satisfy H*P(v) = b');
            PPv = Pv + H \ (b - H*Pv);
            tc.verifyLessThan(norm(PPv - Pv) / norm(Pv), tc.TIGHT_TOL, ...
                'the projection must be idempotent');
            [Q, ~] = qr(full(H)', 0);             % P(v) - v must lie in range(H')
            d = Pv - v;
            tc.verifyLessThan(norm(d - Q*(Q'*d)) / norm(d), 1e-10, ...
                'P(v) - v must be orthogonal to null(H)');
            clearfactors(H);
            tc.verifyEmpty(H.factorcache.ulv, 'clearfactors must drop the stored factors');
            tc.verifyLessThan(norm(H \ b - x) / norm(x), tc.TIGHT_TOL, ...
                'H\\b after clearfactors must refactor and give the same x');
        end

        function test_modified_copy_drops_stored_factors(tc)
            m = 64; n = 128;
            rng(9);
            H = hss(tc.nestedLowRank(m, n, tc.BSIZE, 2), blocksize = tc.BSIZE, tol = 1e-12);
            b = randn(m, 1);
            x = H \ b;
            H2 = H;                               % a copy shares the factors ...
            tc.verifyNotEmpty(H2.factorcache.ulv, 'an unmodified copy may reuse the factors');
            H2.A11.D = 2 * H2.A11.D;              % ... until it is modified
            tc.verifyEmpty(H2.factorcache.ulv, 'modifying H2 must drop its stored factors');
            x2 = H2 \ b;
            xref = pinv(full(H2)) * b;
            tc.verifyLessThan(norm(x2 - xref) / norm(xref), tc.TIGHT_TOL, ...
                'the modified copy must be solved with its own factors');
            tc.verifyEqual(H \ b, x, 'the original must keep its own factors');
        end

        function test_square_solve_stores_factors(tc)
            % square H\b uses the same ULV factorization (no size reduction),
            % keeps its factors, and takes several right-hand sides at once
            n = 64;
            rng(5);
            H = hss(tc.nestedLowRank(n, n, tc.BSIZE, 2), blocksize = tc.BSIZE, tol = 1e-12);
            Hf = full(H);
            B = randn(n, 3);
            X = H \ B;
            tc.verifyLessThan(norm(X - Hf\B) / norm(Hf\B), tc.TIGHT_TOL, ...
                'square H\\B must solve every column');
            tc.verifyNotEmpty(H.factorcache.ulv, 'square H\\b must keep its factors');
            x2 = H \ B(:, 2);                    % one column, from the stored factors
            tc.verifyLessThan(norm(x2 - X(:, 2)) / norm(X(:, 2)), tc.TIGHT_TOL, ...
                'a repeat square solve must reuse the factors');
        end

        function test_minnorm_weighted_leaf_blocks(tc)
            m = 64; n = 128;
            rng(10);
            H = hss(tc.nestedLowRank(m, n, tc.BSIZE, 2), blocksize = tc.BSIZE, tol = 1e-12);
            Hf = full(H);
            b = randn(m, 1);
            w = 10.^(2*(rand(n, 1) - 0.5));
            xref = pinv(Hf * diag(1./w)) * b ./ w;
            tc.verifyLessThan(norm(minnorm(H, b, 'Weight', w) - xref) / norm(xref), tc.TIGHT_TOL, ...
                'minnorm with a diagonal weight must give argmin ||diag(w) x|| s.t. H x = b');
            p = tc.leafcols(H);                   % leaf column block sizes
            Lb = cell(numel(p), 1);
            for i = 1:numel(p)
                Lb{i} = 3*eye(p(i)) + 0.5*randn(p(i));
            end
            Lf = blkdiag(Lb{:});
            xref = Lf \ (pinv(Hf / Lf) * b);
            x = minnorm(H, b, 'Weight', Lb);
            tc.verifyLessThan(norm(x - xref) / norm(xref), tc.TIGHT_TOL, ...
                'minnorm with leaf blocks must give argmin ||L x|| s.t. H x = b');
            tc.verifyEqual(minnorm(H, b, 'Weight', Lb), x, ...
                'a repeat call with the same weight must reuse the stored factors');
        end

        function test_tikhonov_leaf_blocks_any_shape(tc)
            rng(11);
            for sz = [64 128; 64 64; 128 64]'
                m = sz(1); n = sz(2);
                if m > n
                    A = tc.nestedLowRank(n, m, tc.BSIZE, 2)';   % helper builds wide/square only
                else
                    A = tc.nestedLowRank(m, n, tc.BSIZE, 2);
                end
                H = hss(A, blocksize = tc.BSIZE, tol = 1e-12);
                Hf = full(H);
                b = randn(m, 1);
                lam = 0.1;
                p = tc.leafcols(H); q = tc.leafrows(H);
                Lb = cell(numel(p), 1); Sb = cell(numel(q), 1);
                for i = 1:numel(p)
                    Lb{i} = 3*eye(p(i)) + 0.5*randn(p(i));
                    Sb{i} = eye(q(i)) + 0.3*randn(q(i));
                end
                Lf = blkdiag(Lb{:}); Sf = blkdiag(Sb{:});
                xref = [Sf*Hf; lam*Lf] \ [Sf*b; zeros(n, 1)];
                [x, r] = tikhonov(H, b, lam, 'Weight', Lb, 'DataWeight', Sb);
                tc.verifyLessThan(norm(x - xref) / norm(xref), 1e-9, ...
                    sprintf('tikhonov (%dx%d) must minimize ||S(Hx-b)||^2 + lam^2||Lx||^2', m, n));
                tc.verifyLessThan(norm(r - Sf*(Hf*x - b)) / norm(r), 1e-9, ...
                    'tikhonov must return r = S(Hx - b)');
            end
        end

        function test_minnorm_complex(tc)
            m = 64; n = 128;
            rng(6);
            A = tc.nestedLowRank(m, n, tc.BSIZE, 2) + 1i*tc.nestedLowRank(m, n, tc.BSIZE, 2);
            H = hss(A, blocksize = tc.BSIZE, tol = 1e-12);
            Hf = full(H);
            b = randn(m, 1) + 1i*randn(m, 1);
            x = H \ b;
            xref = pinv(Hf) * b;
            tc.verifyLessThan(norm(x - xref) / norm(xref), tc.TIGHT_TOL, ...
                'complex H\\b must be the minimum-norm solution');
        end

        function test_minnorm_backward_stable_illconditioned(tc)
            % 'valid' Gaussian blur, sigma = 2: kappa ~ 1e8. A backward-
            % stable solve gives a backward error ~1e-16 (a root solve by
            % pinv(D)*b gives ~4e-13 here).
            n = 300; sigma = 2; w = 10;
            g = exp(-((-w:w)/sigma).^2/2);
            m = n - 2*w;
            A = zeros(m, n);
            for i = 1:m
                A(i, i:i+2*w) = g;
            end
            H = hss(A, blocksize = tc.BSIZE, tol = 1e-15);
            Hf = full(H);
            rng(7);
            b = randn(m, 1);
            x = H \ b;
            bwd = norm(b - Hf*x) / (norm(Hf)*norm(x) + norm(b));
            tc.verifyLessThan(bwd, 1e-14, ...
                'H\\b must be backward stable on an ill-conditioned matrix');
        end

        function test_minnorm_rank_deficient_errors(tc)
            % The first leaf is 12 x 4 and couples to the rest of H through a
            % rank-2 block, so 10 of its rows live on its 4 columns only: H
            % cannot have full row rank. The solver must stop with a named
            % error instead of returning a vector that does not solve
            % H*x = b. (cutrule places the row cut at 12 and the column cut
            % at 4; the tree has one level, so it is called only at the root.)
            rng(8);
            A = [randn(12, 4), randn(12, 2)*randn(2, 20); ...
                 randn(8, 2)*randn(2, 4), randn(8, 20)];
            H = hss(A, blocksize = 8, tol = 1e-10, cutrule = @(k) 12*(k == 20) + 4*(k == 24));
            tc.verifyEqual(size(H.A11), [12, 4], 'sanity: the first leaf must be 12 x 4');
            tc.verifyError(@() H \ randn(20, 1), 'hss_ulvminnormsolve:rankDeficient', ...
                'rank-deficient H must raise hss_ulvminnormsolve:rankDeficient');
        end
    end


    % ==================================================================
    % 5c. SOLVER ON MATRIX FAMILIES
    %
    % H\b on structured matrix families, built from generators with
    % hss_examples/lib (never through hss()): F1 low rank + block
    % diagonal, F2 random HSS (real, complex, m/n up to 0.95, l+m > n leaves, t = 0 leaves, depth
    % 6, square), F3 Cauchy, F4 convolution, F5 NUDFT (also with a gap on
    % a geometry-aligned tree), F6 Toeplitz, Gaussian blur (kappa ~ 1e8).
    % Every case checks, against the dense minimum-norm solution:
    %   error <= max(1e-11, 100 kappa eps), backward error <= 1e-14,
    %   one stored factor level per tree level (no merging), and that a
    %   repeat solve from the stored factors gives the same x.
    % Plus weighted minnorm / tikhonov with complex data and vector
    % weights, and the improperRanks error.
    % ==================================================================
    methods (Test, TestTags = {'minnorm', 'families'})

        function test_families_generators(tc)
            C = {
                'F1 LR + block diagonal, rank 3',       mp_gen_lrbd(512, 1024, 3, 16, 1)
                'F2 random HSS, rank 10',               mp_gen_random(512, 1024, 10, 32, 2)
                'F2 complex, rank 10',                  mp_gen_random(512, 1024, 10, 32, 3, struct('complex', true))
                'F2 m/n=0.50, rank 16',                 mp_gen_random(512, 1024, 16, 64, 4)
                'F2 m/n=0.90, rank 16 (slack fails)',   mp_gen_random(922, 1024, 16, 64, 4)
                'F2 m/n=0.95, rank 16 (slack fails)',   mp_gen_random(973, 1024, 16, 64, 4)
                'F2 16x24 leaves, rank 12 (l+m > n)',   mp_gen_random(256, 384, 12, 16, 7)
                'F2 rank = leaf rows (t = 0)',          mp_gen_random(256, 512, 16, 16, 5)
                'F2 depth-6 tree, rank 8',              mp_gen_random(1024, 2048, 8, 16, 6)
                'F2 square 1024, rank 10',              mp_gen_random(1024, 1024, 10, 32, 9)
                };
            for c = 1:size(C, 1)
                tc.checkFamilySolve(C{c, 1}, mp_hss_from_generators(C{c, 2}), 20 + c);
            end
        end

        function test_families_kernels(tc)
            K = mp_kernel_cauchy(342, 3); [L, rb, cb] = mp_tree(K.m, K.n, 32);
            tc.checkFamilySolve('F3 interlaced Cauchy', mp_hss_from_generators(mp_hss_kernel(K, L, rb, cb, 1e-12)), 31);
            K = mp_kernel_conv(2048, 2, 'ricker', 3); [L, rb, cb] = mp_tree(K.m, K.n, 32);
            tc.checkFamilySolve('F4 Ricker convolution, stride 2', mp_hss_from_generators(mp_hss_kernel(K, L, rb, cb, 1e-14)), 32);
            rng(11); x5 = sort(((0:767)' + 0.5 + 0.5*(rand(768, 1) - 0.5))/768);
            K = mp_kernel_nudft(x5, 1024); [L, rb, cb] = mp_tree(K.m, K.n, 32);
            tc.checkFamilySolve('F5 NUDFT Cauchy-like (slack fails)', mp_hss_from_generators(mp_hss_kernel(K, L, rb, cb, 1e-12)), 33);
            rng(12); t6 = randn(512+1024-1, 1)/16; tcol = t6(512:-1:1); trow = t6(512:end); tcol(1) = trow(1);
            K = mp_kernel_toeplitz(tcol, trow); [L, rb, cb] = mp_tree(512, 1024, 32);
            tc.checkFamilySolve('F6 Toeplitz -> Cauchy-like', mp_hss_from_generators(mp_hss_kernel(K, L, rb, cb, 1e-12)), 34);
            rng(13); xg = sort(rand(700, 1)); xg = xg(xg < 0.42 | xg > 0.52); N = 1024;
            K = mp_kernel_nudft(xg, N); [L, rb, cb] = mp_tree_aligned(mod(-xg, 1), (0:N-1)'/N, 32);
            tc.checkFamilySolve('F5 with a gap, geometry-aligned tree', mp_hss_from_generators(mp_hss_kernel(K, L, rb, cb, 1e-12)), 35);
        end

        function test_families_gaussian_blur(tc)
            for sig = [1.6 2.0]
                K = mp_kernel_conv(300, 1, 'gauss', sig, ceil(5*sig)); [L, rb, cb] = mp_tree(K.m, K.n, 16);
                tc.checkFamilySolve(sprintf('Gaussian blur sigma=%.1f', sig), ...
                    mp_hss_from_generators(mp_hss_kernel(K, L, rb, cb, 1e-15)), 40 + round(10*sig));
            end
        end

        function test_weighted_tikhonov_complex_and_vectors(tc)
            for cplx = [false true]
                G = mp_gen_random(256, 512, 8, 32, 12, struct('complex', cplx));
                H = mp_hss_from_generators(G);
                A = full(H); [m, n] = size(A);
                rng(49); b = randn(m, 2) + 1i*cplx*randn(m, 2);
                cb = G.cb{end}; rb = G.rb{end}; nl = numel(cb) - 1;
                w = 10.^(2*(rand(n, 1) - 0.5)); s = 1 + rand(m, 1);
                Lb = cell(nl, 1); Sb = cell(nl, 1);
                for i = 1:nl
                    p = cb(i+1) - cb(i); q = rb(i+1) - rb(i);
                    Lb{i} = 3*eye(p) + 0.5*randn(p); Sb{i} = eye(q) + 0.3*randn(q);
                end
                Lf = blkdiag(Lb{:}); Sf = blkdiag(Sb{:}); lam = 0.3;
                tag = sprintf(' (complex=%d)', cplx);
                xr = diag(1./w) * mp_minnorm_dense(A*diag(1./w), b);
                tc.verifyLessThan(mp_rel(minnorm(H, b, 'Weight', w), xr), 1e-11, ['minnorm, vector weight' tag]);
                xr = Lf \ mp_minnorm_dense(A/Lf, b);
                tc.verifyLessThan(mp_rel(minnorm(H, b, 'Weight', Lb), xr), 1e-11, ['minnorm, leaf-block weight' tag]);
                xr = [Sf*A; lam*Lf] \ [Sf*b; zeros(n, 2)];
                [x, r] = tikhonov(H, b, lam, 'Weight', Lb, 'DataWeight', Sb);
                tc.verifyLessThan(mp_rel(x, xr), 1e-10, ['tikhonov, leaf-block L and S' tag]);
                tc.verifyLessThan(mp_rel(r, Sf*(A*x - b)), 1e-10, ['tikhonov residual r = S(Hx - b)' tag]);
                xr = [diag(s)*A; lam*diag(w)] \ [diag(s)*b; zeros(n, 2)];
                tc.verifyLessThan(mp_rel(tikhonov(H, b, lam, 'Weight', w, 'DataWeight', s), xr), 1e-10, ...
                    ['tikhonov, vector L and S' tag]);
            end
        end

        function test_improper_ranks_errors(tc)
            % a leaf row basis of rank 6 on a 4-row leaf: the representation
            % is not proper and the solver must say so
            rng(14);
            G = struct('L', 1, 'rb', {{[0 8], [0 4 8]}}, 'cb', {{[0 24], [0 12 24]}});
            G.D = {randn(4, 12), randn(4, 12)};
            G.U = {{}, {randn(4, 6), randn(4, 6)}};
            G.V = {{}, {orth(randn(12, 6)), orth(randn(12, 6))}};
            G.B12 = {{randn(6)}}; G.B21 = {{randn(6)}};
            H = mp_hss_from_generators(G);
            tc.verifyError(@() H \ randn(8, 1), 'hss_ulvminnormsolve:improperRanks');
        end
    end


    % ==================================================================
    % 6. FULL RECONSTRUCTION  (full(H) -- a second, independent
    %    dense-reconstruction path from extract()/subsref, which every
    %    reconstruction test above goes through via the dense() helper)
    % ==================================================================
    methods (Test, TestTags = {'full'})

        function test_full_matches_dense_multilevel(tc)
            n = tc.N_MED;
            A = tc.cauchy(n);
            H = tc.build(A);
            relerr = norm(A - full(H), 'fro') / norm(A, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'full(H) must match the original matrix');
        end

        function test_full_matches_extract(tc)
            n = tc.N_MED;
            H = tc.build(tc.cauchy(n));
            tc.verifyEqual(full(H), tc.dense(H), 'AbsTol', 1e-10, ...
                'full(H) and H(:,:) (extract) must agree');
        end

        function test_full_small_single_leaf(tc)
            % Regression: full() had no base case for H.isleaf -- a
            % single-leaf HSS matrix (whole matrix <= blocksize, no
            % A11/A12/A21/A22 at all) crashed trying to recurse into
            % nonexistent children.
            n = tc.N_LEAF;
            A = tc.cauchy(n);
            H = hss(A, blocksize = n);
            tc.verifyTrue(logical(H.isleaf), 'sanity: H must be a single leaf here');
            tc.verifyEqual(full(H), A, 'AbsTol', 1e-13, ...
                'full(H) must equal A exactly for a single-leaf HSS matrix');
        end

        function test_full_rectangular(tc)
            m = 32; n = 64;
            A = tc.cauchy(m, n);
            H = hss(A, blocksize = tc.BSIZE);
            relerr = norm(A - full(H), 'fro') / norm(A, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'full(H) must match A for a rectangular HSS matrix');
        end
    end


    % ==================================================================
    % 7. EXTRACT / INDEXING  ( H(rows,cols) -- every reconstruction test
    %    above only ever extracted the WHOLE matrix via dense(); this
    %    exercises genuinely partial ranges )
    % ==================================================================
    methods (Test, TestTags = {'extract'})

        function test_extract_range_within_single_leaf(tc)
            % blocksize=16: rows/cols 3:10 stay entirely inside the leaf
            % spanning indices 1:16 -- no leaf boundary crossed.
            n = tc.N_MED;
            A = tc.cauchy(n);
            H = tc.build(A);
            tc.verifyEqual(H(3:10, 3:10), A(3:10, 3:10), 'AbsTol', 1e-10, ...
                'extraction within a single leaf must match the dense submatrix');
        end

        function test_extract_range_crossing_leaf_boundary(tc)
            % blocksize=16: rows/cols 10:20 straddle the boundary between
            % the two leaves at [1:16] and [17:32].
            n = tc.N_MED;
            A = tc.cauchy(n);
            H = tc.build(A);
            tc.verifyEqual(H(10:20, 10:20), A(10:20, 10:20), 'AbsTol', 1e-10, ...
                'extraction crossing a leaf boundary must match the dense submatrix');
        end

        function test_extract_range_crossing_top_level_split(tc)
            % Rows/cols 30:40 straddle the ROOT's own A11/A22 split (at
            % n/2=32), not just a leaf-pair boundary a level down.
            n = tc.N_MED;
            A = tc.cauchy(n);
            H = tc.build(A);
            tc.verifyEqual(H(30:40, 30:40), A(30:40, 30:40), 'AbsTol', 1e-10, ...
                'extraction crossing the top-level split must match the dense submatrix');
        end

        function test_extract_single_element(tc)
            n = tc.N_MED;
            A = tc.cauchy(n);
            H = tc.build(A);
            tc.verifyEqual(H(5,5), A(5,5), 'AbsTol', 1e-10, ...
                'single-element extraction must match');
        end

        function test_extract_nonrectangular_row_col_ranges(tc)
            % Independently-ranged, differently-sized row/col ranges on a
            % rectangular H.
            m = 32; n = 64;
            A = tc.cauchy(m, n);
            H = hss(A, blocksize = tc.BSIZE);
            tc.verifyEqual(H(5:12, 40:55), A(5:12, 40:55), 'AbsTol', 1e-10, ...
                'independently-ranged row/col extraction must match');
        end

        function test_extract_scrambled_and_duplicate_indices(tc)
            % Regression: extract.m split requested indices across
            % A11/A22 via intersect(), which sorts its result ascending
            % and drops duplicates -- unlike A(idx,idx), which preserves
            % the caller's order and repeats. rowidx/colidx below are
            % deliberately unsorted, contain a repeat, and straddle the
            % top-level split (at 32) so both quadrant-splitting and
            % result-scattering are exercised.
            n = tc.N_MED;
            A = tc.cauchy(n);
            H = tc.build(A);
            rowidx = [45, 10, 45, 3, 30];
            colidx = [5, 50, 20, 50];
            tc.verifyEqual(H(rowidx, colidx), A(rowidx, colidx), 'AbsTol', 1e-10, ...
                'extraction must preserve the caller''s row/col order and duplicates, not sort/dedupe them');
        end

        function test_extract_colon_full_and_partial(tc)
            % Regression: H(:,:) (and H(:,cols)/H(rows,:)) passed the bare
            % ':' char straight through to extract(), which then did char
            % arithmetic on its char code (58) instead of indexing every
            % row/col.
            n = tc.N_MED;
            A = tc.cauchy(n);
            H = tc.build(A);
            tc.verifyEqual(H(:,:), A, 'AbsTol', 1e-10, 'H(:,:) must equal the full dense matrix');
            tc.verifyEqual(H(:, 5:10), A(:, 5:10), 'AbsTol', 1e-10, 'H(:,cols) must equal A(:,cols)');
            tc.verifyEqual(H(5:10, :), A(5:10, :), 'AbsTol', 1e-10, 'H(rows,:) must equal A(rows,:)');
        end
    end


    % ==================================================================
    % 8. TRANSPOSE  ( H.' and H' )
    % ==================================================================
    methods (Test, TestTags = {'transpose'})

        function test_transpose_real_matches_dense(tc)
            n = tc.N_MED;
            A = tc.cauchy(n);
            H = tc.build(A);
            relerr = norm(full(H.') - A.', 'fro') / norm(A, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, 'H.'' must match A.''');
        end

        function test_transpose_rectangular_size(tc)
            m = 32; n = 64;
            H = hss(tc.cauchy(m,n), blocksize = tc.BSIZE);
            Ht = H.';
            tc.verifyEqual(size(Ht), [n, m], 'H.'' must swap dimensions');
        end

        function test_ctranspose_real_matches_transpose(tc)
            % For real data, ' and .' must agree.
            n = tc.N_MED;
            A = tc.cauchy(n);
            H = tc.build(A);
            tc.verifyEqual(full(H'), full(H.'), 'AbsTol', 1e-10, ...
                'for real data, H'' and H.'' must agree');
        end

        function test_ctranspose_complex_conjugates(tc)
            % Regression: with no ctranspose.m, H' on a scalar hss object
            % hit MATLAB's default array-level ctranspose (a no-op for a
            % 1x1 object) and silently returned H completely unchanged,
            % rather than erroring or conjugate-transposing. Now fixed.
            n = tc.N_MED;
            A = tc.cauchy(n) + 1i*tc.cauchy(n);
            H = hss(A, blocksize = tc.BSIZE);
            relerr = norm(full(H') - A', 'fro') / norm(A, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'H'' must match A'' (conjugate transpose) for complex data');
        end

        function test_transpose_complex_does_not_conjugate(tc)
            % H.' must not conjugate any part of the tree (for real data
            % ' and .' agree, so only complex data can catch this).
            n = tc.N_MED;
            A = tc.cauchy(n) + 1i*tc.cauchy(n);
            H = hss(A, blocksize = tc.BSIZE);
            relerr = norm(full(H.') - A.', 'fro') / norm(A, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'H.'' must match A.'' (plain transpose, no conjugation) for complex data');
        end
    end


    % ==================================================================
    % 9. HSS x HSS MULTIPLY  ( H1*H2 -- hss_matmat )
    % ==================================================================
    methods (Test, TestTags = {'matmat'})

        function test_matmat_single_leaf(tc)
            n = tc.N_LEAF;
            A = tc.cauchy(n); B = tc.cauchy(n);
            HA = hss(A, blocksize = n); HB = hss(B, blocksize = n);
            HC = HA*HB;
            tc.verifyEqual(full(HC), A*B, 'AbsTol', 1e-10, ...
                'single-leaf H1*H2 must match A*B');
        end

        function test_matmat_leaf_pair_square(tc)
            n = tc.N_SMALL;   % levelcount=1: root's children are leaves
            A = tc.cauchy(n); B = tc.cauchy(n);
            HA = hss(A, blocksize = tc.BSIZE); HB = hss(B, blocksize = tc.BSIZE);
            HC = HA*HB;
            relerr = norm(full(HC) - A*B, 'fro') / norm(A*B, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'H1*H2 must match A*B at levelcount=1');
        end

        function test_matmat_leaf_pair_rectangular(tc)
            m = 32; k = 32; n = 48;   % share the same partition/levelcount
            A = tc.cauchy(m,k); B = tc.cauchy(k,n);
            HA = hss(A, blocksize = tc.BSIZE); HB = hss(B, blocksize = tc.BSIZE);
            HC = HA*HB;
            relerr = norm(full(HC) - A*B, 'fro') / norm(A*B, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'rectangular H1*H2 must match A*B');
            tc.verifyEqual(size(HC), [m, n], 'H1*H2 size must be [m, n]');
        end

        function test_matmat_leaf_pair_complex(tc)
            n = tc.N_SMALL;
            A = tc.cauchy(n) + 1i*tc.cauchy(n);
            B = tc.cauchy(n) + 1i*tc.cauchy(n);
            HA = hss(A, blocksize = tc.BSIZE); HB = hss(B, blocksize = tc.BSIZE);
            HC = HA*HB;
            relerr = norm(full(HC) - A*B, 'fro') / norm(A*B, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, 'complex H1*H2 must match A*B');
        end

        function test_matmat_deep_tree_square(tc)
            % levelcount>1: exercises the general nested/telescoping-basis
            % recursion (translation matrices above the leaf-pair level),
            % not just the levelcount<=1 direct-compression case.
            n = tc.N_LARGE;
            A = tc.cauchy(n); B = tc.cauchy(n);
            HA = hss(A, blocksize = tc.BSIZE); HB = hss(B, blocksize = tc.BSIZE);
            tc.verifyGreaterThan(HA.levelcount, 1, ...
                'sanity: this case must be deeper than levelcount=1');
            HC = HA*HB;
            relerr = norm(full(HC) - A*B, 'fro') / norm(A*B, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'deep H1*H2 must match A*B');
        end

        function test_matmat_deep_tree_rectangular_complex(tc)
            % Deep, rectangular, complex, and with A's/B's ranks differing
            % at every level (independent random-ish Cauchy data) --
            % stresses the interleaved [A-part,B-part] basis bookkeeping
            % in the off-diagonal construction.
            m = 96; k = 128; p = 80;
            A = tc.cauchy(m,k) + 1i*tc.cauchy(m,k);
            B = tc.cauchy(k,p) + 1i*tc.cauchy(k,p);
            HA = hss(A, blocksize = tc.BSIZE); HB = hss(B, blocksize = tc.BSIZE);
            tc.verifyGreaterThan(HA.levelcount, 1, ...
                'sanity: this case must be deeper than levelcount=1');
            HC = HA*HB;
            relerr = norm(full(HC) - A*B, 'fro') / norm(A*B, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'deep rectangular complex H1*H2 must match A*B');
            tc.verifyEqual(size(HC), [m, p], 'H1*H2 size must be [m, p]');
        end

        % ---- Broad structural sweep ----------------------------------

        function test_matmat_broad_sweep(tc)
            % A grid of shapes/depths/blocksizes/dtypes: every product
            % must both match the dense reference and come back as a
            % genuinely well-formed HSS tree (see verifyHssStructure),
            % not just something that happens to flatten correctly.
            rng(7);
            % Each triple [m k n] multiplies an m-by-k by a k-by-n matrix;
            % all are chosen (checked against hss_constructor's leafLevel
            % rule) so A and B end up with the SAME levelcount at both
            % blocksizes below -- an m/k/n mix that's too lopsided (e.g.
            % [40 40 24]) can legitimately give A and B different depths,
            % which is exactly the mismatch hss_matmat must reject, not
            % something to exercise here.
            shapes = { [16 16 16], [32 32 32], [64 64 64], [128 128 128], ...
                       [48 32 64], [96 128 80], [56 40 72] };
            blocksizes = [8, 16];
            for si = 1:numel(shapes)
                sh = shapes{si};
                m = sh(1); k = sh(2); n = sh(3);
                for bs = blocksizes
                    for cplx = [false true]
                        A = tc.cauchy(m,k); B = tc.cauchy(k,n);
                        if cplx
                            A = A + 1i*tc.cauchy(m,k);
                            B = B + 1i*tc.cauchy(k,n);
                        end
                        HA = hss(A, blocksize = bs);
                        HB = hss(B, blocksize = bs);
                        HC = HA*HB;
                        relerr = norm(full(HC) - A*B, 'fro') / norm(A*B, 'fro');
                        tc.verifyLessThan(relerr, tc.LOOSE_TOL, sprintf( ...
                            'shape [%d %d %d] bs=%d cplx=%d: H1*H2 must match A*B', ...
                            m, k, n, bs, cplx));
                        tc.verifyEqual(size(HC), [m, n], sprintf( ...
                            'shape [%d %d %d] bs=%d cplx=%d: wrong output size', ...
                            m, k, n, bs, cplx));
                        tc.verifyHssStructure(HC);
                    end
                end
            end
        end

        function test_matmat_output_is_valid_hss(tc)
            % The product must be usable exactly like any other hss()
            % output: fed into full(), *, ', and another hss_matmat call.
            n = tc.N_LARGE;
            A = tc.cauchy(n); B = tc.cauchy(n);
            HA = hss(A, blocksize = tc.BSIZE); HB = hss(B, blocksize = tc.BSIZE);
            HC = HA*HB;

            tc.verifyHssStructure(HC);
            tc.verifyEqual(HC.levelcount, HA.levelcount, ...
                'product must share the row-side tree depth');
            tc.verifyTrue(logical(HC.isroot), 'product root must have isroot=true');

            x = randn(n,1);
            relerr = norm(HC*x - (A*B)*x) / norm((A*B)*x);
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'product must support matvec like any other HSS matrix');

            HCt = HC.';
            relerr = norm(full(HCt) - (A*B).', 'fro') / norm(A*B, 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'product must support transpose like any other HSS matrix');

            HD = hss(tc.cauchy(n), blocksize = tc.BSIZE);
            HE = HC*HD;   % product-of-product: HC must be a valid LEFT factor too
            relerr = norm(full(HE) - (A*B)*tc.dense(HD), 'fro') / norm(full(HE), 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, ...
                'a product must itself be usable as an operand of hss_matmat');
        end

        function test_matmat_chaining_associative(tc)
            n = tc.N_MED;
            A = tc.cauchy(n); B = tc.cauchy(n) + n*eye(n); C = tc.cauchy(n);
            HA = hss(A, blocksize = tc.BSIZE);
            HB = hss(B, blocksize = tc.BSIZE);
            HC = hss(C, blocksize = tc.BSIZE);

            left  = (HA*HB)*HC;
            right = HA*(HB*HC);
            dense_ref = A*B*C;

            relerr_left  = norm(full(left)  - dense_ref, 'fro') / norm(dense_ref, 'fro');
            relerr_right = norm(full(right) - dense_ref, 'fro') / norm(dense_ref, 'fro');
            tc.verifyLessThan(relerr_left,  tc.LOOSE_TOL, '(HA*HB)*HC must match A*B*C');
            tc.verifyLessThan(relerr_right, tc.LOOSE_TOL, 'HA*(HB*HC) must match A*B*C');
            tc.verifyHssStructure(left);
            tc.verifyHssStructure(right);
        end

        % ---- H*H' / H*H.'  (Gram matrices -- a primary use case) ------

        function test_matmat_hermitian_gram_square(tc)
            n = tc.N_MED;
            A = tc.cauchy(n) + 1i*tc.cauchy(n);
            H = hss(A, blocksize = tc.BSIZE);
            HG = H*H';
            relerr = norm(full(HG) - A*A', 'fro') / norm(A*A', 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, 'H*H'' must match A*A''');
            herm_gap = norm(full(HG) - full(HG)', 'fro') / norm(full(HG), 'fro');
            tc.verifyLessThan(herm_gap, tc.LOOSE_TOL, 'H*H'' must be (numerically) Hermitian');
            tc.verifyHssStructure(HG);
        end

        function test_matmat_hermitian_gram_tall(tc)
            % The primary use case: forming a Gram matrix from tall data.
            m = 96; n = 32;
            A = tc.cauchy(m,n) + 1i*tc.cauchy(m,n);
            H = hss(A, blocksize = tc.BSIZE);
            HG = H*H';
            relerr = norm(full(HG) - A*A', 'fro') / norm(A*A', 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, 'tall H*H'' must match A*A''');
            tc.verifyEqual(size(HG), [m, m], 'H*H'' of an m x n matrix must be m x m');
            tc.verifyHssStructure(HG);
        end

        function test_matmat_hermitian_gram_wide(tc)
            m = 32; n = 96;
            A = tc.cauchy(m,n) + 1i*tc.cauchy(m,n);
            H = hss(A, blocksize = tc.BSIZE);
            HG = H*H';
            relerr = norm(full(HG) - A*A', 'fro') / norm(A*A', 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, 'wide H*H'' must match A*A''');
            tc.verifyEqual(size(HG), [m, m], 'H*H'' of an m x n matrix must be m x m');
            tc.verifyHssStructure(HG);
        end

        % ---- Extreme / degenerate shapes -------------------------------

        function test_matmat_extreme_aspect_ratios(tc)
            % Aspect ratios far outside the earlier tests: very tall times
            % small, small times very wide, and prime/odd dimensions that
            % never split evenly. Some of these degenerate to a
            % single-leaf HSS matrix (hss_constructor requires BOTH
            % dimensions to clear blocksize before splitting further,
            % independent of hss_matmat) -- that degeneracy is exercised
            % here too, via the existing levelcount==0 base case.
            cases = { ...
                struct('m',200,'k',20,'n',20,'bs',8), ...
                struct('m',20, 'k',20,'n',200,'bs',8), ...
                struct('m',97, 'k',53,'n',71, 'bs',8), ...
                struct('m',300,'k',12,'n',300,'bs',4) ...
            };
            for i = 1:numel(cases)
                c = cases{i};
                A = tc.cauchy(c.m,c.k) + 1i*tc.cauchy(c.m,c.k);
                B = tc.cauchy(c.k,c.n) + 1i*tc.cauchy(c.k,c.n);
                HA = hss(A, blocksize = c.bs); HB = hss(B, blocksize = c.bs);
                HC = HA*HB;
                relerr = norm(full(HC) - A*B, 'fro') / norm(A*B, 'fro');
                tc.verifyLessThan(relerr, tc.LOOSE_TOL, sprintf( ...
                    '[%d %d %d] bs=%d: H1*H2 must match A*B', c.m, c.k, c.n, c.bs));
                tc.verifyEqual(size(HC), [c.m, c.n]);
                tc.verifyHssStructure(HC);
            end
        end

        function test_matmat_gram_from_very_tall_data(tc)
            % The realistic driver for H*H': a Gram/normal-equations matrix
            % formed from very tall, skinny data (aspect ratio 30-60x).
            for mn = {[500 15], [1000 32]}
                m = mn{1}(1); n = mn{1}(2);
                A = tc.cauchy(m,n) + 1i*tc.cauchy(m,n);
                H = hss(A, blocksize = tc.BSIZE);
                HG = H*H';
                relerr = norm(full(HG) - A*A', 'fro') / norm(A*A', 'fro');
                tc.verifyLessThan(relerr, tc.LOOSE_TOL, sprintf( ...
                    '%dx%d: H*H'' must match A*A''', m, n));
                tc.verifyEqual(size(HG), [m, m]);
                tc.verifyHssStructure(HG);
            end
        end

        function test_matmat_transpose_product(tc)
            % Plain (non-conjugate) transpose product H*H.' -- distinct
            % code path from ctranspose (see transpose.m vs ctranspose.m).
            n = tc.N_MED;
            A = tc.cauchy(n) + 1i*tc.cauchy(n);
            H = hss(A, blocksize = tc.BSIZE);
            HG = H*H.';
            relerr = norm(full(HG) - A*A.', 'fro') / norm(A*A.', 'fro');
            tc.verifyLessThan(relerr, tc.LOOSE_TOL, 'H*H.'' must match A*A.'' (no conjugation)');
            tc.verifyHssStructure(HG);
        end

        % ---- Error handling for incompatible operands -----------------

        function test_matmat_size_mismatch_errors(tc)
            HA = hss(tc.cauchy(32,32), blocksize = tc.BSIZE);
            HB = hss(tc.cauchy(48,32), blocksize = tc.BSIZE);   % 48 ~= 32
            tc.verifyError(@() HA*HB, 'hss_matmat:sizeMismatch', ...
                'incompatible inner dimensions must error');
        end

        function test_matmat_levelcount_mismatch_errors(tc)
            HA = hss(tc.cauchy(32,32), blocksize = 16);   % levelcount=1
            HB = hss(tc.cauchy(32,32), blocksize = 8);    % deeper tree
            tc.verifyNotEqual(HA.levelcount, HB.levelcount, ...
                'sanity: the two operands must actually differ in depth');
            tc.verifyError(@() HA*HB, 'hss_matmat:levelcountMismatch', ...
                'mismatched tree depth must error');
        end

        function test_matmat_incompatible_structure_errors(tc)
            % Same overall size AND the same levelcount, but built with
            % different cutrules -- so the shared dimension is split at
            % different points internally. Must be rejected explicitly
            % rather than silently misaligning the recursion (or crashing
            % with an opaque low-level dimension-mismatch error).
            n = 64; bs = 8;
            A = tc.cauchy(n); B = tc.cauchy(n);
            HA = hss(A, blocksize = bs);
            HB = hss(B, blocksize = bs, cutrule = @(k) max(1,floor(k/3)));
            tc.verifyEqual(HA.levelcount, HB.levelcount, ...
                'sanity: both trees must have the same depth');
            tc.verifyNotEqual(size(HA.A11, 2), size(HB.A11, 1), ...
                'sanity: the two operands must actually split differently');
            tc.verifyError(@() HA*HB, 'hss_matmat:incompatibleStructure', ...
                'mismatched split points must error, not silently misalign');
        end
    end


    % ==================================================================
    % 10. REGRESSIONS: size/end, indexing, pieces of the tree, adaptive
    %     ID rank, constructor input checks, rank-deficient root, weight
    %     checks.
    % ==================================================================
    methods (Test, TestTags = {'regressions'})

        function test_size_method_forms(tc)
            A = tc.cauchy(96, 64);
            H = hss(A, blocksize = tc.BSIZE);
            tc.verifyEqual(size(H), [96 64]);
            [m, n] = size(H);
            tc.verifyEqual([m n], [96 64]);
            tc.verifyEqual(size(H, 1), 96);
            tc.verifyEqual(size(H, 2), 64);
            tc.verifyEqual(size(H, 3), 1);
            tc.verifyEqual(length(H), 96);
            tc.verifyEqual(ndims(H), 2);
            tc.verifyFalse(isempty(H));
        end

        function test_end_in_indexing(tc)
            A = tc.cauchy(tc.N_MED, 48);
            H = hss(A, blocksize = tc.BSIZE);
            tc.verifyEqual(H(end, end), A(end, end), 'AbsTol', 1e-10);
            tc.verifyEqual(H(end, 1), A(end, 1), 'AbsTol', 1e-10);
            tc.verifyEqual(H(2:end, end-1), A(2:end, end-1), 'AbsTol', 1e-10);
            tc.verifyEqual(H(end), A(end), 'AbsTol', 1e-10);
        end

        function test_index_out_of_range_errors(tc)
            n = tc.N_MED;
            H = hss(tc.cauchy(n), blocksize = tc.BSIZE);
            tc.verifyError(@() H(n+1, 1), 'hss:index:outOfRange');
            tc.verifyError(@() H(1, n+1), 'hss:index:outOfRange');
            tc.verifyError(@() H(0, 1), 'hss:index:invalid');
            tc.verifyError(@() H(2.5, 1), 'hss:index:invalid');
            tc.verifyError(@() H(-1, 1), 'hss:index:invalid');
            tc.verifyError(@() H(1, 1, 2), 'hss:index:outOfRange');
        end

        function test_index_logical_and_linear(tc)
            n = tc.N_MED;
            A = tc.cauchy(n);
            H = hss(A, blocksize = tc.BSIZE);
            mask = false(n, 1); mask([3 40]) = true;
            tc.verifyEqual(H(mask, 1:3), A(mask, 1:3), 'AbsTol', 1e-10);
            tc.verifyEqual(H(1:5, mask), A(1:5, mask), 'AbsTol', 1e-10);
            tc.verifyEqual(H(5), A(5), 'AbsTol', 1e-10);
            tc.verifyEqual(H([3 70 200]), A([3 70 200]), 'AbsTol', 1e-10);
            tc.verifyEqual(H([3; 70]), A([3; 70]), 'AbsTol', 1e-10);
            tc.verifyEqual(H(:), A(:), 'AbsTol', 1e-10);
            tc.verifyEqual(H(1, :), A(1, :), 'AbsTol', 1e-10);
            tc.verifyEqual(size(H([], 1:3)), [0 3]);
        end

        function test_tree_piece_is_refused(tc)
            n = tc.N_MED;
            A = tc.cauchy(n);
            H = hss(A, blocksize = tc.BSIZE);
            P = H.A22;
            tc.verifyError(@() P(1, 1), 'hss:notRoot');
            tc.verifyError(@() P * ones(size(P, 2), 1), 'hss:notRoot');
            tc.verifyError(@() P \ ones(size(P, 1), 1), 'hss:notRoot');
            % full() of a piece still gives its dense block
            tc.verifyEqual(full(P), A(n/2+1:n, n/2+1:n), 'RelTol', 1e-8);
        end

        function test_id_rank_above_sixty(tc)
            % exact HSS with coupling rank 70: every leaf block row has rank
            % 140, above the old sketch cap (60, or 120 after redraws)
            rng(11);
            A = tc.nestedLowRank(1024, 1024, 256, 70);
            H = hss(A, blocksize = 256, tol = 1e-12);
            x = randn(1024, 1);
            tc.verifyLessThan(norm(H*x - A*x) / norm(A*x), 1e-10);
            tc.verifyEqual(size(H.A11.A12.Z, 2), 140);
        end

        function test_cutrule_must_split_properly(tc)
            tc.verifyError(@() hss(randn(64), blocksize = 8, cutrule = @(k) k-1), 'hss:cutrule');
            tc.verifyError(@() hss(randn(64), blocksize = 8, cutrule = @(k) k/3), 'hss:cutrule');
        end

        function test_sizeA_checked(tc)
            A = tc.cauchy(tc.N_MED);
            tc.verifyError(@() hss(A, blocksize = tc.BSIZE, sizeA = [32 32]), 'hss:sizeA');
            % function-handle input still uses sizeA
            Hf = hss(@(I, J) A(I, J), blocksize = tc.BSIZE, sizeA = size(A));
            tc.verifyEqual(size(Hf), size(A));
            tc.verifyLessThan(norm(full(Hf) - A, 'fro') / norm(A, 'fro'), 1e-10);
            tc.verifyError(@() hss(@(I, J) A(I, J), blocksize = tc.BSIZE), 'hss:sizeA');
        end

        function test_nonfinite_input_errors(tc)
            A = randn(tc.N_MED); A(5, 40) = NaN;
            tc.verifyError(@() hss(A, blocksize = tc.BSIZE), 'hss:nonFinite');
            A(5, 40) = Inf;
            tc.verifyError(@() hss(A, blocksize = tc.BSIZE), 'hss:nonFinite');
        end

        function test_rank_deficient_root(tc)
            % two equal rows in different leaves: the dependence reaches the
            % root. No solution for a generic b -> error; consistent b ->
            % warning and the minimum-norm solution.
            rng(12);
            A = tc.nestedLowRank(256, 512, 32, 5);
            A(200, :) = A(10, :);
            H = hss(A, blocksize = 32);
            tc.verifyError(@() H \ randn(256, 1), 'hss_ulvminnormsolve:inconsistent');
            b = A * randn(512, 1);
            x = tc.verifyWarning(@() H \ b, 'hss_ulvminnormsolve:rankDeficientRoot');
            xr = pinv(A) * b;
            tc.verifyLessThan(norm(x - xr) / norm(xr), 1e-10);
        end

        function test_weight_checks(tc)
            rng(13);
            A = tc.nestedLowRank(256, 512, 32, 5);
            H = hss(A, blocksize = 32);
            nl = numel(tc.leafcols(H));
            L = repmat({eye(64)}, 1, nl); L{3}(5, :) = 0;
            b = randn(256, 1);
            tc.verifyError(@() minnorm(H, b, Weight = L), 'hss:weight:singular');
            w = ones(512, 1); w(7) = NaN;
            tc.verifyError(@() minnorm(H, b, Weight = w), 'hss:weight:nonFinite');
            tc.verifyError(@() tikhonov(H, b, 0.1, DataWeight = [NaN; ones(255, 1)]), 'hss:weight:nonFinite');
        end
    end


    % ==================================================================
    % Shared helper: recursively verify an HSS node's own bookkeeping is
    % internally self-consistent (sizes, index ranges, leaf contents),
    % independent of what numbers it actually stores. Used to confirm
    % hss_matmat's output is a genuinely well-formed HSS tree -- one that
    % every other overloaded operator can walk -- not just something that
    % happens to flatten to the right dense matrix.
    % ==================================================================
    methods
        function checkFamilySolve(tc, name, H, seed)
            % H\b vs the dense minimum-norm solution: error, backward error,
            % one factor level per tree level, and a repeat solve from the
            % stored factors
            A = full(H); [m, ~] = size(A);
            rng(seed); b = randn(m, 1);
            if ~isreal(A), b = b + 1i*randn(m, 1); end
            clearfactors(H);
            x = H \ b;
            xh = H \ b;
            xr = mp_minnorm_dense(A, b);
            sv = svd(A); kappa = sv(1)/sv(end);
            tc.verifyLessThan(mp_rel(x, xr), max(1e-11, 100*kappa*eps), [name ': error vs dense min-norm']);
            tc.verifyLessThan(norm(b - A*x)/(sv(1)*norm(x) + norm(b)), 1e-14, [name ': backward error']);
            tc.verifyEqual(tc.factorDepth(H.factorcache.ulv), H.levelcount, [name ': no tree level merged']);
            tc.verifyEqual(xh, x, [name ': repeat solve from the stored factors']);
        end

        function verifyHssStructure(tc, H)
            % A single-leaf root (whole matrix small enough to be one
            % leaf) leaves Ir/Ic empty -- a pre-existing hss() convention,
            % not something to enforce here.
            if ~isempty(H.Ir)
                tc.verifyEqual(H.Ir(2)-H.Ir(1)+1, size(H, 1), 'Ir must span size(1) rows');
                tc.verifyEqual(H.Ic(2)-H.Ic(1)+1, size(H, 2), 'Ic must span size(2) cols');
            end
            if H.isleaf
                tc.verifyEqual(size(H.D), size(H), 'leaf D size must match node size');
                return
            end
            tc.verifyEqual(size(H.A11, 1)+size(H.A22, 1), size(H, 1), ...
                'child row sizes must sum to parent');
            tc.verifyEqual(size(H.A11, 2)+size(H.A22, 2), size(H, 2), ...
                'child col sizes must sum to parent');
            tc.verifyEqual(H.A11.Ir(2)+1, H.A22.Ir(1), 'child row ranges must be contiguous');
            tc.verifyEqual(H.A11.Ic(2)+1, H.A22.Ic(1), 'child col ranges must be contiguous');
            tc.verifyEqual(size(H.A12), [size(H.A11, 1), size(H.A22, 2)], ...
                'A12 size must be [rows(A11), cols(A22)]');
            tc.verifyEqual(size(H.A21), [size(H.A22, 1), size(H.A11, 2)], ...
                'A21 size must be [rows(A22), cols(A11)]');
            tc.verifyEqual(size(H.A12.Z,2), size(H.A12.lrcomponent,1), ...
                'A12 Z/lrcomponent inner dimension must agree');
            tc.verifyEqual(size(H.A12.Y,1), size(H.A12.lrcomponent,2), ...
                'A12 Y/lrcomponent inner dimension must agree');
            tc.verifyHssStructure(H.A11);
            tc.verifyHssStructure(H.A22);
        end
    end
end