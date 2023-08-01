function H = hss_nudftv(nodes, a, b,N, varargin)
% this code constructs an HSS factorization
% for the Cauchy-like matrix VF^*, where V is the nudft matrix
% and F is a scaled DFT. V is size M by N, where M = length(nodes)
%
% It preserves an interpolative structure so that
% the ULV-like solver in [1] can be used. 
% 
% We use ADI on the disp eqn D(VF^*) - (VF^*)L = a*b, where
% D = diag(nodes), L = diag(w^(2*(1:n))), w = exp(pi*1i/n).
 
block_size = 128; %smallest blocks in col dim.  
tol = hssoption('threshold'); %rel tol. 
if ~isempty(varargin)
    tolcheck = find(strcmp('tol', varargin));
    if ~isempty(tolcheck)
        tol = varargin{tolcheck+1};  
    end
end 

%% Build tree
M = length(nodes); 

nda = angle(nodes); jj = find(nda < 0,1); 
nda(jj:end) = nda(jj:end) + 2*pi; 
rou = 2*pi*(1:N)/N;
% Prepare the tree for the HSS structure -- leaving all the blocks empty:
%H2 = hss_build_hss_tree(M, N, block_size);
H = hss_build_hss_tree_nudft(M, N, nda, rou, block_size);

%build the HSS representation recursively:
H = BuildHSS_iter(H, tol, 0, 0,N, nodes, a, b');
end



%% SUBROUTINES
% 
function H = BuildHSS_iter(H, tol, mh, nh,N,nodes, a, b)

m = size(H, 1);
n = size(H, 2);

if H.leafnode    
    %need to do cols first since col ranks will control ranks. 
    
    % Let's do it for the HSS block row first.
    ridx = (mh+1:mh+m).'; 
    cidx = [(1:nh) (nh+n+1:N)].'; 
    cidxl = circshift(cidx, -nh);
    w = exp(pi*1i/N); 
    % I = interval on unit circle used for determining ADI shift params. 
    %I = [nodes(ridx(1)),nodes(ridx(end)), w.^(2*(cidxl(1))), w.^(2*cidxl(end))];
    I = [w.^(2*(cidxl(1))), w.^(2*cidxl(end))];
    I = [flip(I.*[w.^(-1), w.^(1)]), I];
    [p, ~] = getshifts_adi(I, 'tol', tol); 
    kk = min(length(p), n); %can't go bigger than num cols. 
    % NOTE: this means the block isn't actually low rank!
    % use fadi to build an ID:
    aa = length(ridx);  
    [H.U, J] = vfADI_row(aa, I, ridx, nodes,N, kk, a(ridx)); %N is global size 
    H.ridx = J + mh; %keep track of row indices. 

    %     %DELETE WHEN DONE begin error check: 
%     testrows = buildcauchy(ridx, cidx, N, nodes, a, b);
%     buildcheck = H.U*testrows(J,:);
%     if norm(buildcheck-testrows) > 1e-5
%         error('error line 66 (leaf block row)')
%     end
%     %end error check
   
    % And now do the columns 
    M = length(nodes);
    ridx = [(1:mh) (mh+m+1:M)].'; 
    ridxl = circshift(ridx, -mh); 
    cidx = (nh+1:nh+n).'; 
    %kk = length(J);  
    bb = length(cidx); 
    I = [w.^(2*(cidx(1))), w.^(2*cidx(end))];
    I = [flip(I.*[w.^(-1), w.^(1)]), I];
    [H.V, J] = vfADI_col(bb, I,cidx, N,kk, b(cidx)); %N is global size
    H.cidx = J +nh; %keep track of col indices. 

%     %DELETE WHEN DONE begin error check: 
%     testcols = buildcauchy(ridx, cidx, N, nodes, a, b);
%     buildcheck = testcols(:,J)*(H.V)';
%     if norm(buildcheck-testcols) > 1e-5
%         error('error line 85 (leaf block col)')
%     end
%     %end error check


    
    %H.D = buildcauchy(mh+1:mh+m, nh+1:nh+n, N, nodes, a, b); 
    H.D = vbuildcauchydiags(nodes,mh+1:mh+m, nh+1:nh+n, N);
else
    % Call the constructor recursively on the left and right children.
      H.A11 = BuildHSS_iter(H.A11, tol, mh, nh, N, nodes, a, b);
    
      [mm1, nn] = size(H.A11); 
      [mm2, ~] = size(H.A22);
      mh = mh + mm1; 
      nh = nh + nn; 
      H.A22 = BuildHSS_iter(H.A22,tol, mh,nh,N, nodes, a, b);
      
      %now build glue matrices for level below:
      H.B12 = buildcauchy(H.A11.ridx, H.A22.cidx, N, nodes, a, b); 
      H.B21 = buildcauchy(H.A22.ridx, H.A11.cidx, N, nodes, a, b); 
      
      if H.topnode
        return;
      end
      
      %row translation matrices:
      ridx = [H.A11.ridx; H.A22.ridx]; %subset of rows;
       
      idx1 = nh - nn +1; %col left
      idx2 = nh+nn; %col right
      cidx = [(1:idx1-1) (idx2+1:N)].';
      cidxl = circshift(cidx, -(nh-nn));
      aa = length(ridx); 
      
      w = exp(1i*pi/N); 
      I = [w.^(2*(cidxl(1))), w.^(2*cidxl(end))];
      I = [flip(I.*[w.^(-1), w.^(1)]), I];
      
      [R, J] = vfADI_row(aa, I, ridx, nodes, N,tol, a(ridx));
      aa = length(H.A11.ridx); 
      H.Rl = R(1:aa, :); 
      H.Rr = R(aa+1:end, :); 
      H.ridx = ridx(J); 
      kk = length(J);

    %     %DELETE WHEN DONE begin error check: 
%     testrows = buildcauchy(ridx, cidx, N, nodes, a, b);
%     buildcheck = R*testrows(J,:);
%     if norm(buildcheck-testrows) > 1e-5
%         error('error line 136 (transl. block row)')
%     end
%     %end error check


     
      %now col translation matrices: 
      M = length(nodes);
      idx1 = mh - mm1 +1; %row top
      idx2 = mh+mm2; %row bottom
      ridx = [(1:idx1-1) (idx2+1:M)].';
      cidx = [H.A11.cidx; H.A22.cidx]; 
      bb = length(cidx); 
      %ridxl = circshift(ridx, -(mh-mm)); 
      I = [w.^(2*min(cidx)), w.^(2*max(cidx))];
      I = [flip(I.*[w.^(-1), w.^(1)]), I];
      [R, K] = vfADI_col(bb,I,cidx, N, kk, b(cidx)); 
      H.Wl = R(1:aa, :); 
      H.Wr = R(aa+1:end, :); 
      H.cidx = cidx(K); 


    %     %DELETE WHEN DONE begin error check: 
%     testcols = buildcauchy(ridx, cidx, N, nodes, a, b);
%     buildcheck = testcols(:,K)*R';
%     if norm(buildcheck-testcols) > 1e-5
%         error('error line 162 (transl. block col)')
%     end
%     %end error check
end
    
end


%%
function AA = buildcauchy(J, K, N, nodes, G, L)
%build the cauchy matrix from the given indices:
J = J(:); 
K = K(:); 
w = exp(1i*pi/N); 
a = nodes(J);
b = w.^(2*K); 
A = bsxfun(@minus, a, b.'); 
A = 1./A; 


%apply hadamard product
%[~,r] =size(G); 
k = length(J);
n = length(K); 
%AA = 0; 
%for j = 1:r
AA = spdiags(G(J), 0,k,k)*A*spdiags(conj(L(K)), 0, n,n);
    %AA = AA +B; 
%end

end







