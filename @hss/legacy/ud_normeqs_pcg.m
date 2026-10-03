function x = ud_normeqs_pcg(H,b,options)
%MINNORM_DUAL  Min-norm solution of Hx=b via dual normal equations

arguments
    H
    b
    options.tol = 1e-12
    options.maxits = 2000
    %options.x0 = zeros(H.size(2),1)
end

Ht = H.';

AAT = @(y) H*(Ht*y);

[y,~] = pcg(AAT,b,options.tol,options.maxits);
x = Ht*y;
end