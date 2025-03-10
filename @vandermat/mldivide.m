function X = mldivide(V, B)
%MLDIVIDE Left matrix divide for VANDERMAT.
%   Solves the equation VX = B

% Make sure the typing on V,B are correct
% Flip order if necessary
% Throw errors if type matchup can't be handled
if(~isa(V,'vandermat'))
    X = V \ full(B);
else

    % IMPLEMENT THIS:
    % if( pseudos computed )
    %     use them;
    % end
    
    leftpseudo = construct_pseudoinverses(V);
    X = leftpseudo(B);
end

end
