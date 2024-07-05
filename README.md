# structmats
INSTRUCTIONS FOR SETTING UP PACKAGE WITH HM-TOOLBOX are in ```setup.m```

Computing with Vandermonde, Toeplitz, Hankel and other matrices that have special displacement structures:
This is an ongoing collection of algorithms for computing with matrices that have displacement structures. For now, the repo is limited and construction is ongoing.

A major goal is the implementation of rank-structured direct solvers in MATLAB for solving systems involving highly structured matrices (Toeplitz, Vandermonde, Hankel, etc). Our solvers work by transforming the system into one involving a Cauchy-like matrix, and then capitalizing on the rank structure of the Cauchy-like matrix. For now, solvers for type-II nonuniform discrete Fourier transform problems can be found in the repo heatherw3521/NUDFT. A related solver for square Toeplitz matrices is also included in that repository.  


 




