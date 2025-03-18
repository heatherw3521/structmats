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
% systems with dimensions that are powers of 2 as well as 
% overdetermined systems. 
% 
% Once we have our HSS library installed, this will be fixed.


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
% For rectangular, we use disp struct ZgT-TZ1 = GH', where Zg is the circumshift
% matrix except that Zg(1,end)=g. 
%
% For square, we use disp struct Z1T-TZ_{-1} = GH' (as in Xia's paper). 
 
wn = exp(pi*1i/n); 
wm = exp(pi*1i/m); 

%get Toeplitz generators:
if m==n
    error('will input soon...')
 %if ~(log2(m) == ceil(log2(m))) % check a power of 2
 %    error('structsolv_toeplitz: for now, the fast solver requires dimensions to be powers of 2.')
 %end
 %warning off
 %H = hss('cauchytoeplitz', n, G, L,'tol', tol);
 %warning on
 %x = H\b;
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
 %G = G(:,1); L = L(:,1);
 warning off
 H = hss('nudft', nodes, G, L',modes,'tol', tol, 'toep');
 warning on
%%  
 % nodes = nodes(p);
 % rou = wn.^(2*(1:n)).';
 % C = (1./(nodes-rou.')).*(G*L'); 
 % xx = C\b;

 % test buildcauchy
 %J = (1:m);
 %K = (1:n);
 %AA = buildcauchy(J, K, n, nodes,G, L);

 Ln = urv(H); 
 x = urv_solve(Ln,b);
 x = Qr'*fft(x)/sqrt(n);
 if isreal(tc)&& isreal(tr)
     x = real(x);
 end

 if nargout==1
     varargout = {x};
 else
     varargout = {Ln,p, x};
 end
end

end