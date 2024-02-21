function C = vertcat(A, B)
%HORZCAT Vertical Concatenation for TOEPLITZMAT objects
 
% Concatenating two TOEPLITZMAT objects
if( isa(A, 'toeplitzmat') && isa(B, 'toeplitzmat'))
    
    % Check if sizes are compatible
    if( size(A,2) ~= size(B,2) )
        error( 'TOEPLITZMAT:vertcat:sizemismatch', ...
            'Cannot vertically concatenate Toeplitz matrices of sizes %s -by- %s and %s -by- %s.',...
            ""+size(A,1),""+size(A,2), ""+size(B,1),""+size(B,2) );   
    end
    
    % SPECIAL CASE: Result is also toeplitz
    n = size(B,2);
    m = size(A,1);
    if(m >= n-1 && isequal( flip(A.tc(end-n+2:end)), B.tr(2:end).'))
        C = toeplitzmat([A.tc;B.tc], A.tr);
    else
        aflat = [flip(A.tc(2:end)') A.tr];
        if(isequal( aflat(1:n-1), B.tr(2:end)))
            C = toeplitzmat([A.tc;B.tc], A.tr);
        else
            C = [toeplitz(A.tc, A.tr);toeplitz(B.tc, B.tr)];
        end
    end
    
% Concatenating TOEPLITZMAT with something else
else
    
    if( isempty(A)) % [TOEPLITZMAT [] ]
        C = B;
    elseif ( isempty(B) ) % [ [] TOEPLITZMAT]
        C = A;
    elseif ( isa(A, 'toeplitzmat') ) % Attempt the concatenation 
        C = [toeplitz(A.tc, A.tr);B];
    else
        C = [A;toeplitz(B.tc, B.tr)];
    end
end

end

