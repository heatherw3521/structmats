function h = subsref(T, S)
%SUBSREF gets the part of a Toeplitzmat object by referenced subscript.
%   The result may or may not itself be Toeplitz

if(isequal(S(1).type, '{}'))
    error( "TOEPLITZMAT:subsref:unsupportedtype", ...
        'Subindexing with {} unsupported for TOEPLITZMAT.');
end

% SOME WEIRD STUFF WAS HAPPENING, this was the fix
% TODO: figure out if this is truly needed
if(S(1).type == '.')
    h = builtin('subsref',T,S(1));
    if(size(S,2) >= 2)
        h = subsref(h, S(2:end));
    end
    return;
end

% Get dimensions of T
M = size(T,1); 
N = size(T,2);

mi = S.subs{1}; % Row Indexing
if size(S.subs,2) == 1 % Select rows from the first column
    c = T.tc(mi);
    h = toeplitzmat(c, c(1) ); % This will be a vector
    return;
end
m = size(mi,2);

ni = S.subs{2}; % Column Indexing
n = size(ni,2);

% Simplification cases
consecm = false; consecn = false;
if(isequal(mi,':') && isequal(ni,':'))
    h = T; return;
elseif(isequal(mi,':'))
    mi = 1:M; m = M;
    consecm = true;
elseif(isequal(ni,':'))
    ni = 1:N; n = N;
    consecn = true;
end

if(n==1) % Select rows from the ni(1)-th column
    c = dsample(T, mi, ni(1));
    h = toeplitzmat(c, c(1)); % This will be a vector
    return;
end

% Compute whether the result will be toeplitz
% This takes O(nm) time in the worst case, but that is no worse
% than building an n-by-m matrix!
% Catches a common simpler case in O(n+m) time
consec = @(ki) isequal(ki, ki(1):ki(end));

consecm = consecm || consec(mi); % should short circuit for efficiency!
consecn = consecn || consec(ni);

res_toeps = consecm && consecn; % Ininital check
if(~res_toeps) % Deeper check - if needed
    res_toeps = true;
    for j = -m+2:n-2
        % Get the length of the diagonal then get the diagonal
        if(j <= 0)
            K = min(m+j,n);
            dj = mi(1-j:K-j)-ni(1:K);
        else
            K = min(m,n-j);
            dj = mi(1:K)-ni(1+j:j+K);
        end

        if(nnz(dj == dj(1)) == K)
            continue;
        else
            res_toeps = false;
            break;
        end
    end
end

if(res_toeps)
    h = toeplitzmat( dsample(T,mi, ni(1)) , dsample(T,mi(1),ni) );
else
    h = subsref(full(T), S);
end

end