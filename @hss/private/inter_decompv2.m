function [Z, rows, varargout] = inter_decompv2(A, options)
% one-sided or double-sided interpolative decomposition of matrix A
% A = Z*A(rows,:)                 if orientation = 'rows'
% A = A(:,rows)*Z                 if orientation = 'columns'
% A ≈ Z*A(rows,cols)*Z2           if orientation = 'double-sided'
%
% INPUTS:
%   A           : input matrix (m x n)
%   options.ctype : 'threshold' or 'rank' (default: 'threshold')
%   options.cval  : cutoff threshold or rank value (default: 1e-12)
%   options.orientation : 'rows', 'columns', or 'double-sided'
%
% OUTPUTS:
%   Z, rows     : interpolation matrix and index set
%   varargout   : optional outputs for double-sided case: {cols, Z2}

arguments
    A
    options.ctype = 'threshold'
    options.cval = 1e-12
    options.orientation = 0
end

% determine orientation
if options.orientation
    side = options.orientation;
else
    if size(A,2) >= size(A,1)
        side = 'rows';
    else
        side = 'columns';
    end
end
cutofftype = options.ctype;
cutoffval = options.cval;

% Column ID  (A ≈ A(:,cols)*Z)
if strcmp(side, 'columns')

    % Pivoted QR:  A(:,P) = Q * R
    [~, R, P] = qr(A, 'econ', 'vector');

    % Determine numerical rank or fixed rank
    if strcmp(cutofftype, 'threshold')
        k = rank(A / norm(A), cutoffval);
    else
        k = min(cutoffval, min(size(A)));
    end
    if k == 0
        Z = zeros(0, size(A,2));
        rows = [];
        return
    end

    % Extract R11 and R12
    R11 = R(1:k, 1:k);
    if k < size(R,2)
        R12 = R(1:k, k+1:end);
    else
        R12 = zeros(k, 0);
    end

    % Construct interpolation matrix in permuted ordering
    Z_perm = [eye(k), R11 \ R12];  % size k x n_permuted

    % Undo the permutation: compute inverse permutation
    invP = zeros(1, numel(P));
    invP(P) = 1:numel(P);

    % return Z in original column order
    Z = Z_perm(:, invP);

    % columns selected from A
    rows = P(1:k);                 


% Row ID (A ≈ Z*A(rows,:))
% Implemented by performing column ID on A'
elseif strcmp(side, 'rows')

    [Zt, rows] = inter_decompv2(A', ctype=cutofftype, cval=cutoffval, orientation='columns');
    Z = Zt';


% Double-Sided ID (A ≈ Z*A(rows,cols)*Z2)
elseif strcmp(side, 'double-sided')

    [Z, rows] = inter_decompv2(A, ctype=cutofftype, cval=cutoffval, orientation='rows');
    [Z2, cols] = inter_decompv2(A(rows,:), ctype=cutofftype, cval=cutoffval, orientation='columns');

    varargout{1} = cols;
    varargout{2} = Z2;

else
    error("Orientation must be 'rows', 'columns', or 'double-sided'")
end

end
