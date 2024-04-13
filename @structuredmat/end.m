function r = end(M,k,n)
%END End operator for STRUCTUREDMAT.
%   Gets the size of the k-th dimension of a structuredmat M.
%   This is needed for subindexing with expressions like 'M(2:end)' etc
    r = size(M,k);
end

