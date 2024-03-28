function h = subsref(T, S)
%SUBSREF gets the part of a Toeplitzmat object by referenced subscript.
%   The result may or may not itself be Toeplitz - this is determined by
%   a call to 'structure_retained'

% quick hack: only deal with one-argument inputs in the body
if(size(S,2) > 1)
    s1 = subsref(T, S(1));
    h = subsref(s1, S(2:end));

else

    % Brace indexing unsupported for toeplitzmat objects
    if(isequal(S.type, '{}'))
        error( "TOEPLITZMAT:subsref:unsupportedtype", ...
            'Subindexing with {} unsupported for TOEPLITZMAT.');
    end

    % Getting a property
    if(isequal(S.type,'.'))
        h = builtin('subsref',T,S);

    % Getting a submatrix
    elseif(isequal(S.type,'()'))

        % Extract the row subindexing
        I = S.subs{1}.';
        if(length(S.subs) == 1) % Pulling from the first column
            c = T.tc(I);
            h = toeplitzmat(c, c(1) ); % This will be a vector
            return;
        end

        % Extract the column subindexing
        J = S.subs{2};

        % Parse I,J
        if(isequal(I,':')) 
            I = (1:size(T,1))';
        end

        if(isequal(J,':')) 
            J = 1:size(T,2);
        end

        % Check if the structure is retained
        if(structure_retained(T, I, J))
            h = toeplitzmat( dsample(T,I, J(1)) , dsample(T,I(1),J) );
        else
            h = dsample(T,I,J); % Get the full matrix via subindexing
        end
    
    % Weird attempt at subindexing, (this should never run?)
    else
        error( ["TOEPLITZMAT:subsref:unsupportedtype", ...
            newline, ...
            'Subindexing with %s unsupported for TOEPLITZMAT.'],...
            S.type);
    end
    
end

end