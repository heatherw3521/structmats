function h = diag_apply(f, ftest, fname, C, G)
%DIAG_APPLY Applies functions that preserve circulant structure under
%   certain conditions. It turns out that many of the regular operations,
%   namely the pointwise ones, maintain circulant structure if the input
%   pair is either Circulant-Circulant, Circulant-Toeplitz, or Circulant-scalar

% CASE: Second input is toeplitz, give precedence
if( isequal(class(G), 'toeplitzmat') )
    h = diag_apply(f, ftest, fname, toeplitzmat(C), G);
    return;
end

% CASE: Two circulantmat Objects
if( isa(C, 'circulantmat') && isa(G, 'circulantmat'))
    % Fit check
    if( ~isequal(size(C), size(G)) )
        error( ['CIRCULANTMAT:%s:sizemismatch', newline, ...
                'Cannot %s Circulant matrices of sizes %s -by- %s and %s -by- %s.'],...
                fname, ""+size(C,1),""+size(C,2), ""+size(G,1),""+size(G,2) );       
    end
    h = circulantmat( f(C.tc, G.tc) );
    
% CASE: One Toeplitzmat Object    
elseif ( ftest(G) ) % CIRCULANTMAT f ftype
    if( isscalar(G) )
        h = circulantmat(f(C.tc,G));
    else
        h = f(full(C), G);
    end
elseif ( ftest(C) ) % ftype f CIRCULANTMAT
    if( isscalar(C) )
        h = circulantmat(f(C,G.tc));
    else
        h = f(C, full(G));
    end
    
% ERROR CASE: Unsupported type    
else % Operation with unsupported type
    errid = char("CIRCULANTMAT:"+fname+":unsupportedtype"); 
    error( [errid, newline, ...
            'Undefined function ''%s'' for input arguments of type %s ' ...
            'and %s.'], fname, class(C), class(G));
end

end

