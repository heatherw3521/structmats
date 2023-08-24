# structmats
INSTRUCTIONS FOR SETTING UP PACKAGE WITH HM-TOOLBOX are in ```setup.m```

Computing with Vandermonde, Toeplitz, Hankel and other matrices that have special displacement structures:
For now, the repo is limited. It includes solvers for Toeplitz and type-II nonuniform discrete Fourier transform problems. We will be working on expansions, but please contact me if you have a special request. 

This is an ongoing collection of direct solvers in MATLAB for solving systems involving highly structured matrices (Toeplitz, Vandermonde, Hankel, etc), including a solver for inverse nonuniform discrete Fourier transform (NUDFT) problems. Our solvers work by transforming the system into one involving a Cauchy-like matrix, and then capitalizing on the rank structure of the Cauchy-like matrix. 

The software works with objects constructed in the HM-toolbox (see Massei, Robol, Kressner, 2019). Repo: https://github.com/numpi/hm-toolbox. 




