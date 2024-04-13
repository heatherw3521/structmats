function h = mtimes(C, G)
%MTIMES Matrix Times for CIRCULANTMAT objects
%   This is where the magic happens! Matrix-matrix products can be
%   computed much cheaper when one of the two is a toeplitz matrix.

% CASE: we can use .* if one of the inputs is scalar
if(isscalar(C) || isscalar(G))
    h = C .* G;

% CASE: both inputs non-scalar
elseif( ~isa(C, 'circulantmat') ) % ??? * CIRCULANTMAT
    h = mtimes(G', C')';

else % circulantmat * ???
    % CASE: Bad input sizes
    if(size(C,2) ~= size(G,1) )
        error( ['CIRCULANTMAT:mtimes:sizemismatch', newline, ...
            'Cannot multiply Circulant matrices of sizes %s -by- %s and %s -by- %s.'],...
            ""+size(C,1),""+size(C,2), ""+size(G,1),""+size(G,2) );  
    end

    if(isa(G,'circulantmat'))
        G = full(G);
    end

    d = fft(C.tc);        
    h = ifft(d.*(fft(G)));

end
end

