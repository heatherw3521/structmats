function h = rdivide(T, g)
%./   Pointwise TOEPLITZMAT right divide
%   T ./ g pointwise divides T by g. 
%   T and g can be scalars or TOEPLITZMATs

if ( isa(T, 'double') ) % DOUBLE ./ TOEPLITZMAT
    h = toeplitzmat(T ./ g.tc, T ./ g.tr);

elseif ( isa(g,'double') ) % TOEPLITZMAT ./ DOUBLE
    h = toeplitzmat(T.tc ./ g, T.tr ./ g);

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

