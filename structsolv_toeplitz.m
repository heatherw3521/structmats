%structmat solver: 
function out = structsolv_toeplitz(tc,tr, b, varargin)
%solve a linear system Tx = b, where T is Toeplitz. 
% tc = first col of T
% tr = first row of T
% b = right-hand side
%
% structsolv_toeplitz(tc,tr, b, 'tol', tol) sets accuracy to tol.
% The default setting for tol is 1e-11. 
%
% structsolv_toeplitz(tc, tr, b, 's') uses a transformation that involves a 
% singular displacement equation (see [2]). Otherwise, a nonsingular
% displacement equation is used (see [1]). 
%
% NOTES: Currently, the Toeplitz solver is limited to square systems with
% dimensions that are powers of 2. This will be remedied soon. 

% References: 
% [1] Xia, Xi, Gu. "A Superfast Structured Solver
% for Toeplitz Linear Systems via Randomized Sampling."
% SIMAX, Vol. 33 No 3, p.837-858, 2012.
%
% [2] Wilber, H.D. Ch. 4 in "Computing numerically with rational functions", 
% PhD Dissertation, Cornell Univ., 2021.  

if any(strcmpi(varargin{:}, 's'))
    out = Toeplitz_solve(tc,tr, b, varargin);
else
    out = Toeplitz_solve_ns(tc,tr, b, varargin);
end

end
