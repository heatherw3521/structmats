function varargout = structsolv_toeplitz(tc,tr, b, varargin)
% solve a linear system Tx = b, where T is Toeplitz. 
% tc = first col of T
% tr = first row of T
% b = right-hand side
%
% structsolv_toeplitz(tc,tr, b, 'tol', tol) sets accuracy to tol.
% The default setting for tol is 1e-11. 
%  
% Currently this solver is extremely limited; It can handle square 
% systems with dimensions that are powers of 2 as well as 
% overdetermined (rectangular tall--thin) systems. 
% 
% Once we have our HSS library installed, it will be more versatile.


% References: 
% [1] Xia, Xi, Gu. "A Superfast Structured Solver
% for Toeplitz Linear Systems via Randomized Sampling."
% SIMAX, Vol. 33 No 3, p.837-858, 2012.
%
% [2] Beckermann, Kressner, Wilber. "Compression properties in large 
% Toeptlitz-like matrices", 2025,  https://arxiv.org/abs/2502.09823
% 

%%
% check for tolerance parameter: 
tol = 1e-11; 
if ~isempty(varargin)
    if strcmpi(varargin{1}, 'tol')
        tol = varargin{2}; 
    else
        error('structmats:structsolve_toeplitz:could not parse input')
    end
end

n = length(tr); 
m = length(tc);

%% if the matrix is small, solve directly: 
if n <257 && m < 257
    T = toeplitz(tc,tr); 
    x = T\b;
    varargout = {x};
    return
end

%% Part 1: transform to Cauchy system: 
% For rectangular, we use disp struct ZgT-TZ1 = GH', where Zg is the circumshift
% matrix except that Zg(1,end)=g. 
%
% For square, we use disp struct Z1T-TZ_{-1} = GH' (as in Xia's paper). 
 
wn = exp(pi*1i/n); 
wm = exp(pi*1i/m); 

%get Toeplitz generators:
if m==n
 if ~(log2(m) == ceil(log2(m))) % check a power of 2
     error('structsolv_toeplitz: for now, the fast solver for square systems requires dimensions to be powers of 2.')
 end

%get Toeplitz generators:
 [GG, LL] = toep_gens(tr, tc); 
        
%transform to Cauchy-like generators:
 D0 = spdiags(wn.^((1:n).'-1), 0, n,n);
 G = sqrt(m)*ifft(GG);
 L = sqrt(n)*ifft(D0*LL);
 
 %transform RHS: 
 b = sqrt(m)*ifft(b);

 warning off
 H = hss('cauchytoeplitz', n, G, L,'tol', tol);
 warning on
 %x = H\b;
 Ln = ulv(H);
 x = ulv_solve(Ln, b);
 x = D0'*fft(x)/sqrt(n); 
 if isreal(tr)&& isreal(tc)&& isreal(b)
     x = real(x);
 end
 if nargout==1
     varargout = {x};
 else
     varargout = {Ln, x};
 end

elseif m < n
    error('structmats:structsolv_toeplitz: fast solver for underdetermined systems not yet implemented.')
else
q = exp(1i*2*pi*pi/2/m); % ensures that the nodes are shifted off ROU.
[GG, LL] = toep_gens_rect(tr, tc,q^m); 
        
%transform to Cauchy-like generators:
 Qr = spdiags((wn.^(2*(1:n))).', 0, n,n);
 Ql = spdiags((q.^(1:m)).', 0, m, m);
 G = sqrt(m)*ifft(Ql*GG);
 L = sqrt(n)*ifft(Qr*LL);

 %transform RHS: 
 b = sqrt(m)*ifft(Ql*b); 
 
%% Part 2: HSS 
% build HSS approx to C:
%
 nodes = q*(wm.^(2*(0:m-1))).';
 modes = n;
 warning off
 %permute so that nodes are ordered wrt argument:
 args = mod(angle(nodes), 2*pi);
 [args, p] = sort(args); 
%circshift so nodes are ordered appropriately for subdivisions. 
 kk = find(args > pi/n, 1); 
 p = circshift(p, -kk+1); 
 nodes = nodes(p); % permutation equiv to shifting rows of system
 G = G(p,:);
 b = b(p,:);
 warning off
 H = hss('nudft', nodes, G, L',modes,'tol', tol, 'toep');
 warning off

 Ln = urv(H); 
 x = urv_solve(Ln,b);
 x = Qr'*fft(x)/sqrt(n);
 if isreal(tc)&& isreal(tr) && isreal(b)
     x = real(x);
 end

 if nargout==1
     varargout = {x};
 else
     varargout = {Ln,p, x};
 end
end

end