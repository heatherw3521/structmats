function [pol, zer, r, qpoles, qzeros, q, sigma,sigmavec] = Z3(E,F,m,n, nlawson, damping, zolotol)
%Z Computes degree m numerical solutions to Zolotarev's 3rd and 4th problems given
%   sets E and F.
%
%   The third problem is given by 
%       r = arginf_{r in Rn} sup_{E} |r(z)| / inf_{F} |r(z)|,
%   where Rn is the space of type (n,n) rational functions. Loosely stated,
%   we want to make r small on E while keeping it big on F.
%
%   The fourth problem is given by 
%       q = arginf_{q in Rn} ||q-sign_{E/F}||,
%   where Rn is the space of type (n,n) rational functions. Concretely, we
%   want the rational function of type (n,n) which best uniformly
%   approximates the sign_{E/F}.
%   note: sign_{E/F} = 1 on E, -1 on F
%
%   ARGUMENTS:
%   E : Set on which r is small max
%   F : Set on which r has large min
%   m : Degree of rational function
%
%   THIS CODE IS TAKEN DIRECTLY FROM
%   "Computation of Zolotarev rational functions"
%   By Nick Trefethen and Heather Wilber
%   https://arxiv.org/pdf/2408.14092
arguments
    E,F,m

    n = m
    nlawson = 200;
    damping = 0.95;
    zolotol = -1;
end

[q,qpoles,qzeros,zj,fj,wj,~,~,~,tauvec] = Z4(E,F,m,n,true,nlawson,damping,zolotol);
r = @(z) (1+q(z))./(1-q(z)); % Approximate transform, works for tau << 1
[~,~,zer] = prz(zj,fj+1,wj);
[~,~,pol] = prz(zj,fj-1,wj);
% pol = pol(r(pol)>1e12);
% zer = zer(r(zer)<1e-12);
sigma = max(abs(r(E)))/min(abs(r(F)));
sigmavec = tauvec .^2; % Approximately true for tau << 1

end

