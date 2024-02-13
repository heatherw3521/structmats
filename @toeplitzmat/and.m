function h = and(T,g)
%&  Pointwise logical and for TOEPLITZMAT objects
%   T and g can be TOEPLITZMATs or scalars
%   behaves the same as & for matrices in MATLAB

% Check that T,g are actually logical values

classg = class(g);
if(isa(g,'toeplitzmat'))
   classg = class(g.tc);
end
if(~islogical(T) || ~islogical(g))
    error( 'TOEPLITZMAT:and:invalidinputtype', ...
        ['Undefined function ''logical and'' for input arguments of type %s and %s'], ...
        class(T.tc), classg);
end

ftest = @(x) (islogical(x)&&~isa(x,'toeplitzmat'));
h = toepcompare(T, g, @(x,y) x & y, "&",ftest);
end