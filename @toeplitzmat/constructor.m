function T = constructor(T, varargin)
%CONSTRUCTOR The main TOEPLITZMAT constructor
%   A Toeplitz matrix is determined by it's first row and first column.
%   As such, the input will be of the forms:
%
%   CASE 1: |varargin| = 1
%       (a): The input is a CIRCULANTMAT. We cast it to a TOEPLITZMAT.
%       (b): The input is a single vector. To be consistent with MATLAB, we
%            interpret this vector to be both the first column and the
%            first row -- i.e. we create a symmetric TOEPLITZMAT with first
%            row equal to the input vector
%
%   Case 2: varargin is of size 2, containing two vectors.
%           In this case, we are given the first column and the first row.
%           To be consistent with MATLAB, we assume the first vector is the
%           column and the second is the row. The corresponding TOEPLITZMAT
%           is then constructed, after checking that the first-elements of
%           the vectors agree.
%
%   ERROR CASES: There are a few different ways for a user to provide
%               invalid input.
%           (i): |varargin| = 0
%                It does not quite make sense to have an empty constructor
%                for TOEPLITZMAT.
%          (ii): |varargin| > 2
%                We do not support tensors, and so it does not make sense
%                to handle more than 2 inputs for a TOEPLITZMAT
%         (iii): First element mismatch
%                The first elements of the first row and first column of
%                any matrix need to match (this is the top-left element).
%                Thus, we cannot make sense of two input vectors with
%                differing first elements.
%          (iv): Bad input type
%                If the input type does not fit into one of the valid CASES
%                outlined above, then we have no way to make sense of it.
%                For instance, if a string is passed in, we cannot turn it
%                into a TOEPLITZMAT in a meaningful way.
%
%   Here is an example of constructing a toeplitzmat:
%
%   col = [1 2 3];
%   row = [1 4 5];
%   T = toeplitzmat(col,row);
%
%   The code above creates a toeplitzmat representing
%       [1 4 5]
%   T = [2 1 4]
%       [3 2 1]
%   Of course, only the first row and column are stored.

% ERROR (i): Too few args
if(numel(varargin) == 0)
    error(['STRUCTMATS:TOEPLITZMAT:constructor:toofewargs', newline, ...
       'Too few (0) arguments passed in to construct Toeplitz matrix']);  
    
% ERROR (ii)
elseif( numel(varargin) > 2)
    error(['STRUCTMATS:TOEPLITZMAT:constructor:toomanyargs', ...
            newline, ...
            'Too many (%s) arguments passed in to construct Toeplitz matrix'],...
            ""+numel(varargin));
end


% CASE 1: |varargin| = 1
if( numel(varargin) == 1)
    v = varargin{1};

    % (a): Input is a circulantmat, cast it to TOEPLITZMAT
    if( isa(v,'circulantmat') )
        T = toeplitzmat(varargin{1}.tc, varargin{1}.tr);
    
    % (b): Input is a numeric vector, create a symmetric TOEPLITZMAT
    elseif( isvector(v) && (isnumeric(v) || islogical(v)))
        T = toeplitzmat(varargin{1},varargin{1});

    % ERROR (iv): Bad input type
    else
        error(['STRUCTMATS:TOEPLITZMAT:constructor:badinput', ...
            newline, ...
            'Cannot create a toeplitzmat out of single input of type %s'],...
            class(v));
    end

% ERROR (ii): Too many inputs
elseif( numel(varargin) > 2)
    error(['STRUCTMATS:TOEPLITZMAT:constructor:toomanyargs', ...
            newline, ...
            'Too many (%s) arguments passed in to construct Toeplitz matrix'],...
            ""+numel(varargin));

% CASE 2: |varargin| = 2
else
    % Get the first row and column
    tc = varargin{1};
    tr = varargin{2};
    
    % ERROR (iv): Bad input types
    if( ~isvector(tc) || ~isvector(tr))
        error(['STRUCTMATS:TOEPLITZMAT:constructor:badinput', ...
                newline, ...
                'At least one of the inputs was not a vector']);        
    elseif( (~isnumeric(tc) && ~islogical(tc)) || ~isa(tr, class(tc)))
        error(['STRUCTMATS:TOEPLITZMAT:constructor:inputtypemismatch', ...
                newline, ...
                'Inputs pair of types %s and %s is not supported for TOEPLITZMAT.'],...
                class(tc), class(tr));
    end
        
    % Flip-flop tc,tr to ensure that they are a column and a row, resp.
    % i.e. tc should be m-by-1 and tr 1-by-n
    tc = reshape(tc,[],1);
    tr = reshape(tr,1,[]);
    
    % ERROR (iii): Inconsistent top-left element
    if(tc(1) ~= tr(1))
        error(['STRUCTMATS:TOEPLITZMAT:constructor:badinput', ...
                newline, ...
                'Given first column and row have different first element.']);
    end
    
    % All is good - finish construction
    T.tc = tc;
    T.tr = tr;

    % Get Sylvester matrix equation quantities
    m = length(tc); n = length(tr);

    % Form A and B, the circumshift matrices of the right size
    % Here, we make the choice of putting '1' in the top right entry
    A = speye(m);
    A = [A(end,:);A(1:end-1,:)];
    B = speye(n);
    B = [B(end,:);B(1:end-1,:)];
    T.A = A;
    T.B = B;

    % Compute the right hand side
    % Only the first row and last column are nonzero
    y  = [flip(tc).' tr(2:end)];
    dr = y(1:n) - [tr(2:end) tr(1)];
    dc = [y(end-m+1) y(end:-1:end-m+2)].' - tc;
    u1 = [1;zeros(m-1,1)];
    v1 = dr.';
    u2 = [0;dc(2:end)];
    v2 = [zeros(n-1,1);1];
    
    % Form rank-2 right hand side
    T.U = [u1 u2];
    T.V = [v1 v2];
end
end

