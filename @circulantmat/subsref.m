function h = subsref(C, S)
%SUBSREF gets the part of a Toeplitzmat object by referenced subscript.
%   The result may or may not itself be Toeplitz

% This will allow all the code below to just consider 1-element input!
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
    h =  subsref@toeplitzmat(C, S); % USE SUPER
end

end

