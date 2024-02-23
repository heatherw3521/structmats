function h = mtimes(T, G)
%MTIMES Matrix Times for Toeplitzmat objects
%   This is where the magic happens! Matrix-matrix products can be
%   computed much cheaper when one of the two is a toeplitz matrix.

% CASE: we can use .* if one of the inputs is scalar
if(isscalar(T) || isscalar(G))
    h = T .* G;

% CASE: both inputs non-scalar
elseif( ~isa(T, 'toeplitzmat') ) % ??? * Toeplitzmat
    h = mtimes(G', T')';

else % toeplitzmat * ???
    % CASE: Bad input sizes
    if(size(T,2) ~= size(G,1) )
        error( 'TOEPLITZMAT:mtimes:sizemismatch', ...
            'Cannot multiply Toeplitz matrices of sizes %s -by- %s and %s -by- %s.',...
            ""+size(T,1),""+size(T,2), ""+size(G,1),""+size(G,2) );  
    end
    N = size(T,2); 
    M = size(T,1);

    B = toeplitz(G);    
    % RECURSE! (to deal with rectangular matrices
    if (N == 1)
        h = T.tr * toeplitz(G);
    elseif (M ==1)
        h = T.tc * toeplitz(G);
    elseif(N == M)
        h = toep_times(T, toeplitz(G) );
    elseif N > M
        h = toep_times( T(:,1:M), B(1:M,:) ) + ...
            T(:,M+1:end) * B(M+1:end,:);
    else
        h = [toep_times(T(1:N,:), B) ; ...
            T(N+1:end,:) * B];
    end
end
end

