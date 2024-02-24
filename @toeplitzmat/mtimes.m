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

    if(isa(G,'toeplitzmat'))
        G = toeplitz(G);
    end

    CT = circulantmat([T.tc;0;flip(T.tr(2:end)).']);
    y = CT * [G;zeros(M,size(G,2))];
    h = y(1:M,:);
end
end

