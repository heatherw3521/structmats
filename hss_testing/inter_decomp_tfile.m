function [Z,rows,varargout] = inter_decomp_tfile(A,options)
% one sided interpolative decompositon of matrix A (mxn)
% A = Z*A(rows,:) if rows
% A = A(:,cols)*Z if columns

arguments
    A;
    options.ctype = 'threshold';
    options.cval = 1e-12;
    options.orientation = 0;
end


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

% if we want a column ID do row ID on A'
if strcmp(side,'columns')
    % Interpolative decomp of matrix A such that the decomp has rank k
    % uses QR directly on A
    % QR
    [~,R,P] = qr(A,'econ','vector');

    if strcmp(cutofftype,'threshold')
        k = rank(A/norm(A),cutoffval);
    else
        k = min(cutoffval,min(size(A,1),size(A,2)));
    end
    R_k = R(1:k,1:k);
    cols = P(1:k);
    if k<4
        disp(cond(R_k.' * R_k))
    end
    Z = (R_k.' * R_k)\(A(:,cols)'*A);
    rows = cols;
elseif strcmp(side,'rows')
    [Zt, rows] = inter_decomp_tfile(A',ctype = cutofftype,cval = cutoffval,orientation = 'columns');
    Z = Zt';
elseif strcmp(side,'double-sided')
    [Z, rows] = inter_decomp_tfile(A,ctype = cutofftype,cval = cutoffval,orientation ='rows');
    [Z2,cols] = inter_decomp_tfile(A(rows,:),ctype = cutofftype,cval = cutoffval,orientation ='columns');
    % should output Z, rows, cols, Z2 such that A = Z*A(rows,cols)*Z2
    varargout{1} = cols;
    varargout{2} = Z2;
else
    error("Orientation must be 'rows','columns' or 'double-sided'")
end

end


