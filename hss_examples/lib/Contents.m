% HSS_EXAMPLES/LIB  Shared helpers for testing and experimenting with @hss.
%
% Generators -> hss objects (never form the dense matrix)
%   hss_from_generators  - assemble an @hss object from generators
%   hss_to_generators    - read the generators out of an @hss object
%   cluster_tree         - cluster boundaries of the tree hss() builds
%   cluster_tree_aligned - geometry-aligned tree for nonuniform points
%
% Matrix families (generator form or entry evaluators)
%   gen_lowrank_blockdiag - F1: low rank + block diagonal
%   gen_random            - F2: exact random HSS, rank k on every node
%   kernel_cauchy         - F3: interlaced Cauchy kernel
%   kernel_conv           - F4: decimated 1-D convolution (deblurring)
%   conv_apply            - exact y = A*x for a kernel_conv operator
%   kernel_nudft          - F5: Cauchy-like form of the type-II NUDFT
%   kernel_toeplitz       - F6: Cauchy-like form of a Toeplitz matrix
%   hss_from_kernel       - nested-ID HSS construction from an entry evaluator (O(n) proxy mode)
%
% Generator edits
%   gen_scale      - diag(s_rows) * H * diag(s_cols)
%   gen_rightblock - H * blkdiag(M_1, ..., M_q)
%   gen_augment    - Tikhonov augmentation [H, lam*I]
%
% Decompositions
%   id_rows, id_cols - deterministic interpolative decompositions
%   pivqr            - economy column-pivoted QR
%
% Dense references and checks
%   minnorm_dense          - backward-stable dense minimum-norm solution
%   weighted_minnorm_solve - weighted minimum norm via @hss/minnorm
%   tikhonov_solve         - Tikhonov via @hss/tikhonov
%   manufactured_rhs       - right-hand side with a known minimum-norm solution
%   nullspace_fraction     - null-space fraction ||P_null x|| / ||x||
%   relative_error         - relative error
%
% Timing
%   time_median - median wall time of a function
%   time_solve  - first vs repeat H\b timings
%
% Baselines (matvec only)
%   cgne_minnorm         - min-norm solution by CG on the dual normal equations
%   trust_region_minnorm - trust-region / projected-gradient iteration
