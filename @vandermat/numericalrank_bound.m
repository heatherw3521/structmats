function r = numericalrank_bound(V, eps)
%NUMERICALRANK_BOUND Gets a bound on the epsilon-rank of a Vandermonde
%   matrix V.
arguments
    V;
    eps = 1e-12; % Default tolerance
end

% Analytical bound derived in https://arxiv.org/abs/1609.09494 by 
% Beckerman and Townsend (2019)
r = 2*(1+ceil((4*log(8*floor(V.n/2)/pi)*log(4/eps)) / (pi^2)) );
% This seems like a worst-case bound that can be much improved by solving
% the discrete problem

end

