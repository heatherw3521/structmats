function h = plus(T, g)
%+   Plus for TOEPLITZMAT objects
%  F + G adds F and G. F and G can be scalars, TOEPLITZMAT, or matrices

if ( ~isa(T, 'toeplitzmat') ) % ??? + TOEPLITZMAT
    h = plus(g, T);

elseif ( isempty(g) ) % TOEPLITZMAT + []
    h = T;

elseif ( isa(g,'double') ) % TOEPLITZMAT + DOUBLE
    h = toeplitzmat(T.tc + g, T.tr + g);

elseif ( ~isa(g,'toeplitzmat') ) % TOEPLITZMAT + ???
    error( 'TOEPLITZMAT:plus:unknown', ...
        ['Undefined function ''plus'' for input arguments of type %s ' ...
        'and %s.'], class(T), class(g));
 
else                           % TOEPLITZMAT + TOEPLITZMAT
    %Size check
    if( ~isequal(size(T), size(g)) )
        error( 'TOEPLITZMAT:plus:sizemismatch', ...
            'Cannot add Toeplitz matrices of sizes %s -by- %s and %s -by- %s.',...
            ""+size(T.tc,1),""+size(T.tc,2), ""+size(g.tc,1),""+size(g.tr,2) );       
    end
    
    h = toeplitzmat(T.tc + g.tc, T.tr + g.tr);
end

