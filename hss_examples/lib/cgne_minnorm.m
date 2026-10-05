function x = cgne_minnorm(H,b,options)
%CGNE_MINNORM  Minimum-norm solution of H*x = b (H wide, full row rank) by
%   conjugate gradients on the dual normal equations H*H.'*y = b, x = H.'*y.
%   Matrix-vector products only; a baseline for the direct solver H\b.

arguments
    H
    b
    options.tol = 1e-12
    options.maxits = 2000
end

Ht = H.';

AAT = @(y) H*(Ht*y);

[y,~] = pcg(AAT,b,options.tol,options.maxits);
x = Ht*y;
end