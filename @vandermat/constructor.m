function V = constructor(V,varargin)
%CONSTRUCTOR The main VANDERMAT constructor
%   A Vandermonde matrix is determined by it's second column and total
%   number of columns.
%   As such, the input will be of the forms:
%
%   CASE 1: |varargin| = 1
%       The input is a single vector. To be consistent with MATLAB, we
%       interpret that the resulting Vandermonde matrix should be square
%       and further that the given vector is the second column of the
%       matrix.
%
%   Case 2: varargin is of size 2, containing a vector and a positive integer.
%       In this case, we are given the second column and number of columns.
%       The corresponding VANDERMAT is constructed, warning the user
%       if it is singular..
%
%   ERROR CASES: There are a few different ways for a user to provide
%               invalid input.
%           (i): |varargin| = 0
%                It does not quite make sense to have an empty constructor
%                for VANDERMAT.
%          (ii): |varargin| > 2
%                A Vandermonde matrix is resolved by the two pieces of
%                information we are stored (2nd column + # of columns), so
%                there is nothing to do with 3+ inputs.
%          (iii): Bad input types
%                If the input type does not fit into one of the valid CASES
%                outlined above, then we have no way to make sense of it.
%                For instance, if a string is passed in, we cannot turn it
%                into a VANDERMAT in a meaningful way.
%
%   Here is an example of constructing a toeplitzmat:
%
%   x = [1 2 3 4 5];
%   n = 4;
%   V = vandermat(x, n);
%
%   The code above creates a vandermat representing
%       [1 1 1  1  ]
%   V = [1 2 4  16 ]
%       [1 3 9  27 ]
%       [1 4 16 64 ]
%       [1 5 25 125]
%   For efficiency, we are storing only x and n.

% ERROR (i): Too few args
if(numel(varargin) == 0)
    error(['STRUCTMATS:VANDERMAT:constructor:toofewargs', newline, ...
       'Too few (0) arguments passed in to construct Vandermonde matrix']);  
    
% ERROR (ii)
elseif( numel(varargin) > 2)
    error(['STRUCTMATS:VANDERMAT:constructor:toomanyargs', ...
            newline, ...
            'Too many (%s) arguments passed in to construct Toeplitz matrix'],...
            ""+numel(varargin));
end

% CASE 1: |varargin| = 1
if( numel(varargin) == 1)
    x = varargin{1};

    % (a): Input is a numeric vector, create a square VANDERMAT
    if( isvector(x) && (isnumeric(x) || islogical(x)))
        if(size(x,2) > 1) % let's store x as a column vector
            x = x.';
        end
        V.x = x;
        V.n = size(x,1);
    % ERROR (iii): Bad input type
    else
        error(['STRUCTMATS:VANDERMAT:constructor:badinput', ...
            newline, ...
            'Cannot create a vandermat out of single input of type %s'],...
            class(x));
    end

% ERROR (ii): Too many inputs
elseif( numel(varargin) > 2)
    error(['STRUCTMATS:VANDERMAT:constructor:toomanyargs', ...
            newline, ...
            'Too many (%s) arguments passed in to construct Vandermonde matrix'],...
            ""+numel(varargin));

% CASE 2: |varargin| = 2
else
    % Get the second column and number of columns
    x = varargin{1};
    n = varargin{2};
    
    % ERROR (iv): Bad input types
    if( ~isvector(x) )
        error(['STRUCTMATS:VANDERMAT:constructor:badinput', ...
                newline, ...
                'The first input was not a vector.']);
    elseif(n < 1 || mod(n,1)~=0 )
        error(['STRUCTMATS:VANDERMAT:constructor:badinput', ...
                newline, ...
                'The second input must be a positive integer.']);
    end
        
    % Transpose x to make sure it's a column vector.
    x = reshape(x,[],1);
    
    % All is good - finish construction
    V.x = x;
    V.n = n;

    % Get the displacement structure matrices
    V.A = spdiags(x,0,numel(x),numel(x));
    B = speye(n);
    V.B = [B(end,:);B(1:end-1,:)];
    V.u = x.^n-1;
    V.v = [zeros(n-1,1);1];
end

end

