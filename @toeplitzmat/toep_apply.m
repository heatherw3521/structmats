function h = toep_apply(f, ftest, fname, T, G)
%TOEP_APPLY Applies functions that preserve toeplitz structure under
%   certain conditions. It turns out that many of the regular operations,
%   namely the pointwise ones, maintain Toeplitz structure if the input
%   pair is either Toeplitz-Toeplitz or Toeplitz-scalar

% CASE: Two Toeplitzmat Objects
if( isa(T, 'toeplitzmat') && isa(G, 'toeplitzmat'))
    %Size check
    if( ~isequal(size(T), size(G)) )
        error( 'TOEPLITZMAT:%s:sizemismatch', fname, ...
            'Cannot %s Toeplitz matrices of sizes %s -by- %s and %s -by- %s.',...
            fname, ""+size(T.tc,1),""+size(T.tc,2), ""+size(G.tc,1),""+size(G.tr,2) );       
    end
    h = toeplitzmat( f(T.tc, G.tc), f(T.tr, G.tr) );
    
% CASE: One Toeplitzmat Object    
elseif ( ftest(G) ) % TOEPLITZMAT f ftype
    if( isscalar(G) )
        h = toeplitzmat(f(T.tc,G), f(T.tc,G));
    else
        h = f(toeplitz(T.tc,T.tr), G);
    end
elseif ( ftest(T) ) % ftype f TOEPLITZMAT
    if( isscalar(T) )
        h = toeplitzmat(f(T,G.tc), f(T,G.tr));
    else
        h = f(T, toeplitz(G.tc,G.tr));
    end
    
% ERROR CASE: Unsupported type    
else % Operation with unsupported type
    error( "TOEPLITZMAT:"+fname+":unknown",...
        ['Undefined function ''%s'' for input arguments of type %s ' ...
        'and %s.'], fname, class(T), class(G));
end

end

