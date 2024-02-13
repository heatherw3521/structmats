function h = toepcompare(T, g, f, fname)
%TOEPCOMPARE Pointwise comparison for TOEPLITZ objects
%   T and g can be scalars or TOEPLITZMATs
%   f is a binary operation that maps to 0,1. 
%   fname is f's name for error reporting purposes.
%   Really, this script is just for code reuse with lt,gt,le,ge,neq,eq

if ( isa(T, 'double') ) % ??? f TOEPLITZMAT
    h = toeplitzmat(f(T,g.tc), f(T,g.tr)); %rephrase the question

elseif ( isa(g,'double') ) % TOEPLITZMAT f DOUBLE
    h = toeplitzmat(f(T.tc, g), f(T.tr, g));

elseif ( ~isa(g,'toeplitzmat') || ~isa(T,'toeplitzmat')) % TOEPLITZMAT f ???
    error( 'TOEPLITZMAT:%s:unknown', fname, ...
        ['Undefined function ''%s'' for input arguments of type %s ' ...
        'and %s.'], fname, class(T), class(g));
 
else                           % TOEPLITZMAT f TOEPLITZMAT
    %Size check
    if( ~isequal(size(T), size(g)) )
        error( 'TOEPLITZMAT:%s:sizemismatch', fname, ...
            'Cannot compare (%s) Toeplitz matrices of sizes %s -by- %s and %s -by- %s.',...
            fname, ""+size(T.tc,1),""+size(T.tc,2), ""+size(g.tc,1),""+size(g.tr,2) );       
    end
    
    h = toeplitzmat( f(T.tc, g.tc), f(T.tr, g.tr) );
end

