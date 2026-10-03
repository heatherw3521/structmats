function [x,tau_hist] = ud_tr_projgrad(H,b,x0,options)
% tr method
%   x0 must satisfy H*x0 = b

arguments
    H 
    b 
    x0 
    options.maxit = 30
    options.inner_max = 30
    options.alpha = 1
    options.tol = 1e-12
end


tau = norm(x0);
x   = x0;

tau_hist = tau;

for k = 1:options.maxit

    % ---- solve inner trust-region problem ----
    x = trust_region_inner(H,b,tau,x,options.inner_max);

    r = H*x - b;
    g = norm(H.'*r);   % rho'(tau) = -||H.'*r||

    if g < options.tol
        return
    end

    % ---- gradient step in tau ----
    tau = tau - g/tau;
    tau_hist(end+1) = tau;
end
end


function x = trust_region_inner(H,b,tau,x,iters)

rhs = H.'*b;

% Unconstrained candidate
x_unc = cgsolve(@(z) H.'*(H*z), rhs, x);

if norm(x_unc) <= tau
    x = x_unc;
    return
end

% Otherwise solve (H.'H + lambda I)x = H.'b
lambda_lo = 0;
lambda_hi = 1;

while true
    x = cgsolve(@(z) H.'*(H*z) + lambda_hi*z, rhs, x);
    if norm(x) <= tau
        break
    end
    lambda_hi = 2*lambda_hi;
end

for j = 1:iters
    lambda = 0.5*(lambda_lo + lambda_hi);
    x = cgsolve(@(z) H.'*(H*z) + lambda*z, rhs, x);

    if norm(x) > tau
        lambda_lo = lambda;
    else
        lambda_hi = lambda;
    end
end
end

function x = cgsolve(Afun,b,x0)
if nargin < 3, x0 = zeros(size(b)); end
[x,~] = pcg(Afun,b,1e-8,200,[],[],x0);
end
