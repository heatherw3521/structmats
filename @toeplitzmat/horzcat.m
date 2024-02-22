function C = horzcat(varargin)
%HORZCAT Horizontal Concatenation for TOEPLITZMAT objects
 
if(size(varargin,2) == 1)
    C = varargin{1};
elseif(size(varargin,2) > 2)
    C = [varargin{1} horzcat(varargin{2:end})]; 
else
    A = varargin{1};
    B = varargin{2};

    % Concatenating two TOEPLITZMAT objects
    if( isa(A, 'toeplitzmat') && isa(B, 'toeplitzmat'))

        % Check if sizes are compatible
        if( size(A,1) ~= size(B,1) )
            error( 'TOEPLITZMAT:horzcat:sizemismatch', ...
                'Cannot horizontally concatenate Toeplitz matrices of sizes %s -by- %s and %s -by- %s.',...
                ""+size(A,1),""+size(A,2), ""+size(B,1),""+size(B,2) );   
        end

        % SPECIAL CASE: Result is also toeplitz
        m = size(B,1);
        n = size(A,2);
        if(n >= m-1 && isequal( flip(A.tr(end-m+2:end)), B.tc(2:end).'))
            C = toeplitzmat(A.tc, [A.tr B.tr]);
        else
            aflat = [flip(A.tc(2:end)') A.tr];
            if(isequal( aflat(end-m+2:end), B.tc(2:end).'))
                C = toeplitzmat(A.tc, [A.tr B.tr]);
            else
                C = [toeplitz(A.tc, A.tr) toeplitz(B.tc, B.tr)];
            end
        end

    % Concatenating TOEPLITZMAT with something else
    else

        if( isempty(A)) % [TOEPLITZMAT [] ]
            C = B;
        elseif ( isempty(B) ) % [ [] TOEPLITZMAT]
            C = A;
        elseif ( isa(A, 'toeplitzmat') ) % Attempt the concatenation 
            C = [toeplitz(A.tc, A.tr) B];
        else
            C = [A toeplitz(B.tc, B.tr)];
        end
    end
end
end

