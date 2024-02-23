function h = subsref(C, S)
%SUBSREF gets the part of a Toeplitzmat object by referenced subscript.
%   The result may or may not itself be Toeplitz

% This will allow all the code below to consider 1-element input!
if(size(S,2) >= 2)
    h = subsref(subsref(C,S(1)), S(2:end));
end

% SOME WEIRD STUFF WAS HAPPENING, this was the fix
% TODO: figure out if this is truly needed
if(isequal(S.type, '.'))

    if(isequal(S.subs, 'tr'))
        h = C.tc.';
    else
        h =  subsref@toeplitzmat(C, S); % USE SUPER
    end
elseif(isequal(S.type, '()'))
    
    mi = S.subs{1};
    if(size(S.subs,2) == 1 || ~isequal(mi, S.subs{2}))
        h =  subsref@toeplitzmat(C, S); % USE SUPER
    else
        h = circulantmat(C.tc(1:size(mi,2)));
    end

end

end

