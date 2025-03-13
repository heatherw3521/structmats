function [q,qpoles,qzeros,zj,fj,wj,tau,numerator,denominator,tauvec] = Z4(E,F,m,n, forZ3, nlawson, damping, zolotol)
%Z4 Computational solution to Zolotarev's sign approximation problem
arguments
    E,F,m,n

    forZ3 = false;
    nlawson=200;
    damping = 0.95;
    zolotol = -1;
end

fEF = [-ones(size(E)); ones(size(F))];
[q,qpoles,qzeros,zj,numerator,denominator,fj,wj,tauvec] = RationalApprox([E;F],fEF,m,n,true,forZ3,nlawson, damping,1e-13, zolotol);
err = fEF-q([E;F]); err = err(~isnan(err));
tau = norm(err,inf);
end

