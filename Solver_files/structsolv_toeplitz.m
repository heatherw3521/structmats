function varargout = structsolv_toeplitz(tc,tr, b, varargin)
% solve a linear system Tx = b, where T is Toeplitz. 
% tc = first col of T
% tr = first row of T
% b = right-hand side
%
% structsolv_toeplitz(tc,tr, b, 'tol', tol) sets accuracy to tol.
% The default setting for tol is 1e-11. 
%  
% Currently this solver is extremely limited; can handle square 
% systems with dimensions that are powers of 2. Once we have our
% HSS library installed, this will be fixed.

% References: 
% [1] Xia, Xi, Gu. "A Superfast Structured Solver
% for Toeplitz Linear Systems via Randomized Sampling."
% SIMAX, Vol. 33 No 3, p.837-858, 2012.
%
% [2] Beckermann, Kressner, Wilber. "Compression properties in large 
% Toeptlitz-like matrices" https://arxiv.org/abs/2502.09823
% 

%%
% check for tolerance parameter: 
tol = 1e-11; 
if ~isempty(varargin)
    if strcmpi(varargin{1}, 'tol')
        tol = varargin{2}; 
    else
        error('could not parse input')
    end
end

n = length(tr); 
m = length(tc);
block_size = 256;

%% if the matrix is small, solve directly: 
if n <257 && m < 257
    T = toeplitz(tc,tr); 
    x = T\b;
    varargout = {x};
    return
end

%% Part 1: transform to Cauchy system: 

N = (1:n).'; % nodes are w.^(2*(M-1))
M = (1:m).'; 
w1 = exp(pi*1i/n); 
w2 = exp(pi*1i/m); 

%get Toeplitz generators:
[GG, LL] = toep_gens_rect(tr, tc); 
        
%transform to Cauchy-like generators:
 D0 = spdiags(w1.^(-(N-1)), 0, n,n);
 if m==n
    G = sqrt(m)*ifft(GG);
    L = sqrt(n)*ifft(D0'*LL);
 else
    G = sqrt(m)*ifft(GG)*w1;
    L = sqrt(n)*ifft(D0'*LL);
 end

 %transform RHS: 
 b = sqrt(m)*ifft(b); 
 
%% Part 2: HSS 
% build HSS approx to C:
%
if m==n
 if ~(log2(m) == ceil(log2(m))) % check a power of 2
     error('structsolv_toeplitz: for now, the fast solver requires dimensions to be powers of 2.')
 end
 warning off
 H = hss('cauchytoeplitz', n, G, L,'tol', tol);
 warning on
 x = H\b;
else
 error('structsolv_toeplitz: for now, the fast solver works on square systems with dimensions that are powers of 2')
 nodes = (w2.^(2*(M-1)))*w1;
 modes = n;
 warning off
 %permute so that nodes are ordered wrt argument:
 args = mod(angle(nodes), 2*pi);
 [args, p] = sort(args); 
%circshift so nodes are ordered appropriately for subdivisions. 
 kk = find(args > pi/n, 1); 
 p = circshift(p, -kk+1); 
 nodes = nodes(p); % D = diag(nodes)
 G = G(p,:);
 b = b(p,:);

 warning off
 H = hss('nudft', nodes, G, L',modes,'tol', tol, 'toep');
 warning on
%%
   % aa = w2.^(2*(0:m-1))*w1; 
   % bb = w1.^(2*N).';
   % C = bsxfun(@minus, nodes, bb); 
   % C = 1./C; 
   % H2 = (G*L').*C;
%%
%x2 = H\b; 
 Ln = urv(H); 
 x = urv_solve(Ln,b);
end
%% Part 4:
% transform x back: 
 
 x = D0*fft(x)/sqrt(n); 
 if nargout==1
     varargout = {x};
 else
     varargout = {Ln,pp, x};
 end

end