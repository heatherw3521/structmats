function h = diag_apply(f, ftest, fname, H, G)
%DIAG_APPLY Applies functions that preserve hankel structure under
%   certain conditions. It turns out that many of the regular operations,
%   namely the pointwise ones, maintain Toeplitz structure if the input
%   pair is either Toeplitz-Toeplitz or Toeplitz-scalar

% CASE: Two Hankelmat Objects
if( isa(H, 'hankelmat') && isa(G, 'hankelmat'))
    % Fit check
    if( ~isequal(size(H), size(G)) )
        error( 'HANKELMAT:%s:sizemismatch', fname, ...
            'Cannot %s Toeplitz matrices of sizes %s -by- %s and %s -by- %s.',...
            fname, ""+size(H,1),""+size(H,2), ""+size(G,1),""+size(G,2) );       
    end
    h = hankelmat( f(H.hc, G.hc), f(H.hr, G.hr) );
    
% CASE: One Hankelmat Object    
elseif ( ftest(G) ) % HANKELMAT f ftype
    if( isscalar(G) )
        h = hankelmat(f(H.hc,G), f(H.hr,G));
    else
        h = f(hankel(H), G);
    end
elseif ( ftest(H) ) % ftype f HANKELMAT
    if( isscalar(H) )
        h = hankelmat(f(H,G.hc), f(H,G.hr));
    else
        h = f(H, hankel(G));
    end
    
% ERROR CASE: Unsupported type    
else % Operation with unsupported type
    errid = char("HANKELMAT:"+fname+":unsupportedtype"); 
    error( errid, ...
        ['Undefined function ''%s'' for input arguments of type %s ' ...
        'and %s.'], fname, class(H), class(G));
end

end

