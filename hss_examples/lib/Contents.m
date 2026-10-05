% HSS_EXAMPLES/LIB  Shared helpers for testing and experimenting with @hss.
%
% Generators -> hss objects (never form the dense matrix)
%   mp_hss_from_generators - assemble an @hss object from generators
%   mp_hss_to_generators   - read the generators out of an @hss object
%   mp_tree                - cluster boundaries of the tree hss() builds
%   mp_tree_aligned        - geometry-aligned tree for nonuniform points
%
% Matrix families (generator form or entry evaluators)
%   mp_gen_lrbd        - F1: low rank + block diagonal
%   mp_gen_random      - F2: exact random HSS, rank k on every node
%   mp_kernel_cauchy   - F3: interlaced Cauchy kernel
%   mp_kernel_conv     - F4: decimated 1-D convolution (deblurring)
%   mp_conv_apply      - exact y = A*x for an mp_kernel_conv operator
%   mp_kernel_nudft    - F5: Cauchy-like form of the type-II NUDFT
%   mp_kernel_toeplitz - F6: Cauchy-like form of a Toeplitz matrix
%   mp_hss_kernel      - nested-ID HSS construction from an entry evaluator (O(n) proxy mode)
%
% Generator edits
%   mp_gen_scale      - diag(s_rows) * H * diag(s_cols)
%   mp_gen_rightblock - H * blkdiag(M_1, ..., M_q)
%   mp_gen_augment    - Tikhonov augmentation [H, lam*I]
%
% Decompositions
%   mp_id_rows, mp_id_cols - deterministic interpolative decompositions
%   mp_pivqr               - economy column-pivoted QR
%
% Dense references and checks
%   mp_minnorm_dense    - backward-stable dense minimum-norm solution
%   mp_weighted_minnorm - weighted minimum norm via @hss/minnorm
%   mp_tikhonov         - Tikhonov via @hss/tikhonov
%   mp_manufactured     - right-hand side with a known minimum-norm solution
%   mp_nullfrac         - null-space fraction ||P_null x|| / ||x||
%   mp_rel              - relative error
%
% Timing
%   mp_timeit     - median wall time of a function
%   mp_time_solve - first vs repeat H\b timings
%
% Baselines (matvec only)
%   ud_normeqs_pcg - min-norm solution by CG on the dual normal equations
%   ud_tr_projgrad - trust-region / projected-gradient iteration
