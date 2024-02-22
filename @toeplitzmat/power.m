function h = power(T, g)
%.^  Pointwise power for TOEPLITZMAT objects
%  T .^ g raies T to the g, pointwise. 
%  T and g can be scalars or TOEPLITZMATs

if ( isa(T, 'double') ) % DOUBLE .^ TOEPLITZMAT
    if(isscalar(T)) %exponentiating to a number
        h = toeplitzmat(g.tc .^ T, g.tr .^ T);
    else %exponentiating to a matrix
        h = toeplitz(g.tc,g.tr) .^ T;
    end

elseif ( isa(g,'double') ) % TOEPLITZMAT .^ DOUBLE
    if(isscalar(g)) %exponentiating to a number
        h = toeplitzmat(T.tc .^ g, T.tr .^ g);
    else %exponentiating to a matrix
        h = toeplitz(T.tc,T.tr) .^ g;
    end

% TOEPLITZMAT ./ ??? or ??? ./ TOEPLITZMAT
elseif ( ~isa(g,'toeplitzmat') || ~isa(T,'toeplitzmat'))
    error( 'TOEPLITZMAT:power:unknown', ...
        ['Undefined function ''power'' for input arguments of type %s ' ...
        'and %s.'], class(T), class(g));
 
else                           % TOEPLITZMAT .^ TOEPLITZMAT
    %Size check
    if( ~isequal(size(T), size(g)) )
        error( 'TOEPLITZMAT:power:sizemismatch', ...
            'Cannot .^ matrices of sizes %s -by- %s and %s -by- %s.',...
            ""+size(T.tc,1),""+size(T.tc,2), ""+size(g.tc,1),""+size(g.tr,2) );       
    end
    
    h = toeplitzmat(T.tc .^ g.tc, T.tr .^ g.tr);
end

