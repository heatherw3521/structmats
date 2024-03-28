function T = constructor(T, varargin)
%CONSTRUCTOR The main TOEPLITZMAT constructor
%   A Toeplitz matrix is determined by it's first row and first column.
%   As such, the input will be of the forms:
%
%   CASE 1: |varargin| = 1
%       (a): The input is a CIRCULANT. We cast it to a TOEPLITZMAT.
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
                numel(varargin));

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
        if(size(tc,2) > 1)
            tc = tc.';
        end
        if(size(tr, 1) > 1)
            tr = tr.';
        end
        
        % ERROR (iii): Inconsistent top-left element
        if(tc(1) ~= tr(1))
            error(['STRUCTMATS:TOEPLITZMAT:constructor:badinput', ...
                    newline, ...
                    'Given first column and row have different first element.']);
        end
        
        % All is good - finish construction
        T.tc = tc;
        T.tr = tr;
    end
end

