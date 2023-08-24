%structmat solver: 
function varargout = structsolv_toeplitz(tc,tr, b, varargin)
% solve a linear system Tx = b, where T is Toeplitz. 
% tc = first col of T
% tr = first row of T
% b = right-hand side
%
% structsolv_toeplitz(tc,tr, b, 'tol', tol) sets accuracy to tol.
% The default setting for tol is 1e-11. 
% 
%
% NOTES: Currently, the Toeplitz solver is limited to square systems. 
% This will be remedied soon. 

% References: 
% [1] Xia, Xi, Gu. "A Superfast Structured Solver
% for Toeplitz Linear Systems via Randomized Sampling."
% SIMAX, Vol. 33 No 3, p.837-858, 2012.
%
% [2] Wilber, H.D. Ch. 4 in "Computing numerically with rational functions", 
% PhD Dissertation, Cornell Univ., 2021.  


%for now restrict to square: 
if ~(length(tc)==length(tr)) 
    error("error: for now, we only support square Toeplitz matrices.")
end

if any(strcmpi('s',varargin)) %this is not yet ready. Will fail
    out = Toeps_sing(tc,tr, b, varargin{:});
else %this should be fine
    if nargout ==1
        out = Toeps(tc,tr, b, varargin{:});
        varargout = {out};
    else
        [L, x] = Toeps(tc,tr, b, varargin{:});
        varargout = {L,x};
    end
end

end
