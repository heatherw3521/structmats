function c = wrapper_INUDFT(x,N,b,opts)
% WRAPPER_INUDFT  fast direct solver (Cauchy+HSS) 1D NUDFT type 2 lin sys
%
% c = wrapper_INUDFT(x,N,b) solves Ac=b where A is the 1D type 2 (forward) NUDFT
%  matrix with nodes x and N modes. Returns c a length-N column vector of
%  complex Fourier coeffs. NU pts x are on 2pi-periodic domain.
%
%  Converts to Vandermonde form (nodes on unit circle) then calls Epperly+Wilber
%   fast direct solver.
%
% c = wrapper_INUDFT(x,N,b,opts) controls various options such as:
%   optstol - solver requested relative tolerance
%
% For test: see test_iNUDFT

% Barnett 11/18/22
if nargin<4, opts=[]; end
if ~isfield(opts,'tol'), opts.tol = 1e-6; end

% convert to Vandermonde (nodes on unit circle; mode freqs 0:N-1 not -N/2:N/2-1)
ucnodes = exp(1i*x(:));
Noffset = floor(N/2);
rhs = exp(1i*Noffset*x(:)).*b(:);     % prephase the RHS (checked at N=2^10)
try
  c = INUDFT(ucnodes,N,rhs,'tol',opts.tol);        % fast solver, trap errors
catch ME
  fprintf('\tINUDFT failed!\n\t%s\n',ME.message);
  c = nan(N,1);
end
