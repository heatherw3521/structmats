function H = constructor(H, varargin)
%CONSTRUCTOR The main HANKELMAT constructor
%   A Hankel matrix is most conveniently determined by it's first column
%   and last row. As such, the input will be of the forms:
%
%   CASE 1: |varargin| = 1, containing a vector
%           The input is a single vector. To be consistent with MATLAB, we
%           assume that the last row is all zeros, except for the
%           bottom-left element (which is determined by the first column). 
%           We also assume the matrix to be constructed is square.
% 
%   Case 2: varargin is of size 2, containing two vectors.
%           In this case, we are given the first column and the last row.
%           To be consistent with MATLAB, we assume the first vector is the
%           column and the second is the row. The corresponding TOEPLITZMAT
%           is then constructed, after checking that the last element of
%           the first column agrees with the first element of the last row.
%
%   ERROR CASES: There are a few different ways for a user to provide
%               invalid input.
%           (i): |varargin| = 0
%                It does not quite make sense to have an empty constructor
%                for HANKELMAT.
%          (ii): |varargin| > 2
%                We cannot construct a HANKELMAT with more than two inputs
%         (iii): Bottom left element mismatch
%                The first element of the first column and last element of
%                the last row need to match for any matrix (this is the
%                bottom left element). Thus, we cannot make sense of two 
%                input vectors with differing first-last element pair.
%          (iv): Bad input type
%                If the input type does not fit into one of the valid CASES
%                outlined above, then we have no way to make sense of it.
%                For instance, if a string is passed in, we cannot turn it
%                into a HANKELMAT in a meaningful way.
%
%   Here is an example of constructing a hankelmat:
%
%   firstcol = [1 2 3];
%   lastrow = [3 4 5];
%   H = toeplitzmat(col,row);
%
%   The code above creates a toeplitzmat representing
%       [1 2 3]
%   H = [2 3 4]
%       [3 4 5]
%   Of course, only the last row and first column are stored.

%ERROR (i): Too few args
if(numel(varargin) == 0)
    error('STRUCTMATS:HANKELMAT:constructor:toofewargs', newline, ...
            'Too few (0) arguments passed in to construct Hankel matrix');  
        
% ERROR (ii)
elseif( numel(varargin) > 2)
    error(['STRUCTMATS:HANKELMAT:constructor:toomanyargs', ...
            newline, ...
            'Too many (%s) arguments passed in to construct Hankel matrix'],...
            ""+numel(varargin));
end

% Get the first column
hc = varargin{1};

% ERROR (iv): Bad input type (non-vector)
if(~isvector(hc))
    error('STRUCTMATS:HANKELMAT:constructor:badinput', newline, ...
            'Input for first-column is not a vector.');

% ERROR (iv): Bad input type (uninterpretable)
elseif((~isnumeric(hc) && ~islogical(hc)))
    error(['STRUCTMATS:HANKELMAT:constructor:badinput', newline, ...
            'Cannot create a Hankel matrix with input first-column ', ...
            'of type %s.'], class(hc));
end

% transpose hc if needed
if(size(hc, 2) > 1)
    hc = hc.';
end

% CASE 1: single vector input
if( numel(varargin) == 1)
    H.hc = hc;
    H.hr = [hc(end) zeros(1,size(hc,1)-1)];

% CASE 2: two vectors input
else
hr = varargin{2}; % Get the last row

% ERROR (iv): Bad input type (non-vector)
if(~isvector(hr))
    error('STRUCTMATS:HANKELMAT:constructor:badinput', newline, ...
            'Input for last-row is not a vector.');

    % ERROR (iv): Bad input type (mismatch)
elseif( ~isa(hr, class(hc)))
    error(['STRUCTMATS:HANKELMAT:constructor:badinput', newline, ...
            'Inputs pair of types %s and %s is not supported for HANKELMAT.'],...
            class(hc), class(hr));    

% ERROR (iii): Bottom-left element mismatch
elseif(hc(end) ~= hr(1))
    error('STRUCTMATS:HANKELMAT:constructor:firstelementmismatch', ...
            'Given first column and row have different bottom left element.');
end

% Transpose hr if needed
if(size(hr,1) > 1)
    hr = hr.';
end

% All is good - finish construction
H.hc = hc;
H.hr = hr;

end
end

