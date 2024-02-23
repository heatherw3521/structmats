function h = subsref(T,S)
%SUBSREF gets the part of a Toeplitzmat object by referenced subscript.
%   The result may or may not itself be Toeplitz

mi = S.subs{1}; % Row Indexing
if size(S.subs,2) == 1
    h = toeplitzmat(mi, T.tr);
    return;
end

m = size(mi,2);
ni = S.subs{2}; % Column Indexing
n = size(ni,2);

if(isequal(mi, ':') && isequal(ni, ':'))
    h = T;
elseif(isequal(mi,':')
    
elseif(isequal(ni,:))
    h = subsref(toeplitz(T), S);
end


%% TO DO: Make this function smart enough to detect when the result
%   will be Toeplitz

end

