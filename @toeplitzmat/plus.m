function h = plus(T, g)
%+   Plus for TOEPLITZMAT objects
%  T + g adds T and g. 
%  T and g can be scalars or TOEPLITZMATs

if ( ~isa(T, 'toeplitzmat') ) % ??? + TOEPLITZMAT
    h = plus(g, T);

elseif ( isempty(g) ) % TOEPLITZMAT + []
    h = T;

elseif ( isa(g,'double') ) % TOEPLITZMAT + DOUBLE
    if(isscalar(g)) %adding a number
        h = toeplitzmat(T.tc + g, T.tr + g);
    else %adding a matrix
        h = g+toeplitz(T.tc,T.tr);
    end
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

