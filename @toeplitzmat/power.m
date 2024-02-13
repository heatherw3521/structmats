function h = power(T, g)
%.^  Pointwise power for TOEPLITZMAT objects
%  T .^ g raies T to the g, pointwise. 
%  T and g can be scalars or TOEPLITZMATs

if ( ~isa(T, 'toeplitzmat') ) % ??? .^ TOEPLITZMAT
    h = power(g, T);

elseif ( isa(g,'double') ) % TOEPLITZMAT .^ DOUBLE
    h = toeplitzmat(T.tc .^ g, T.tr .^ g);

elseif ( ~isa(g,'toeplitzmat') ) % TOEPLITZMAT .^ ???
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

