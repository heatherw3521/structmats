function h = rdivide(T, g)
%./   Pointwise TOEPLITZMAT right divide
%   T ./ g pointwise divides T by g. 
%   T and g can be scalars or TOEPLITZMATs

if ( isa(T, 'double') ) % DOUBLE ./ TOEPLITZMAT
    if(isscalar(T)) %dividng by a number
        h = toeplitzmat(g.tc ./ T, g.tr ./ T);
    else %dividing by a matrix
        h = toeplitz(g.tc,g.tr) .^ T;
    end

elseif ( isa(g,'double') ) % TOEPLITZMAT ./ DOUBLE
    if(isscalar(g)) %dividing by a number
        h = toeplitzmat(T.tc ./ g, T.tr ./ g);
    else %dividing by a matrix
        h = toeplitz(T.tc,T.tr) ./ g;
    end

% TOEPLITZMAT ./ ??? or ??? ./ TOEPLITZMAT
elseif ( ~isa(g,'toeplitzmat') || ~isa(T,'toeplitzmat') )
    error( 'TOEPLITZMAT:rdivide:unknown', ...
        ['Undefined function ''rdivide'' for input arguments of type %s ' ...
        'and %s.'], class(T), class(g));
 
else                       % TOEPLITZMAT ./ TOEPLITZMAT
    %Size check
    if( ~isequal(size(T), size(g)) )
        error( 'TOEPLITZMAT:rdivide:sizemismatch', ...
            'Cannot elementwise-divide Toeplitz matrices of sizes %s -by- %s and %s -by- %s.',...
            ""+size(T.tc,1),""+size(T.tc,2), ""+size(g.tc,1),""+size(g.tr,2) );       
    end
    
    h = toeplitzmat(T.tc ./ g.tc, T.tr ./ g.tr);
end

