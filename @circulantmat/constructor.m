function C = constructor(C, varargin)
%CONSTRUCTOR The main CIRCULANTMAT constructor
%   A Circulant matrix is determined by it's first column.
%   As such, the input will be of the forms:
%
%   CASE 1: |varargin| = 1, containing a vector
%           The input is the first column of the circulant matrix -- just
%           store it. This determines the entire matrix, but currently, we
%           also want to store the first row. This is due to MATLAB
%           handling subsref in an odd way.
%
%   ERROR CASES: There are a few different ways for a user to provide
%               invalid input.
%           (i): |varargin| = 0
%                It does not quite make sense to have an empty constructor
%                for CIRCULANTMAT.
%          (ii): |varargin| > 1
%                A circulant matrix only needs 1 input: the first column.
%         (iii): Bad input type
%                If the input type does not fit into one of the valid CASES
%                outlined above, then we have no way to make sense of it.
%                For instance, if a string is passed in, we cannot turn it
%                into a CIRCULANTMAT in a meaningful way.
%
%   Here is an example of constructing a circulantmat:
%
%   v = [1 2 3];
%   C = toeplitzmat(v);
%
%   The code above creates a circulantmat representing
%       [1 3 2]
%   C = [2 1 3]
%       [3 2 1]
%   Of course, only the first column and row are stored.
%   (this will be improvewd in the future so that only the first column is
%   stored)

% ERROR (i)
if( numel(varargin) == 0)
    error('STRUCTMATS:CIRCULANTMAT:constructor:toofewargs', ...
            newline, ...
            'Too few (0) arguments passed in to construct Circulant matrix');

% ERROR (ii)
elseif( numel(varargin) >= 2)
    error(['STRUCTMATS:CIRCULANTMAT:constructor:toomanyargs', ...
            newline, ...
            'Too many (%s) arguments passed in to construct Circulant matrix'],...
            numel(varargin));
end

% Get the first column
tc = varargin{1};

% ERROR (iii)
if(~isvector(tc))
    error('STRUCTMATS:CIRCULANTMAT:constructor:badinput', ...
            'Input is not a vector.');
elseif(~islogical(tc) && ~isnumeric(tc))
    error(['STRUCTMATS:CIRCULANTMAT:constructor:badinput', ...
            'Cannot create a CIRCULANTMAT from input of type %s'], ...
            class(tc));
end

% Flip-flop tc to ensure that it is m-by-1
if(size(tc,2) > 1)
    tc = tc.';
end

% All is good - finish construction
C.tc = tc;
C.tr = [tc(1) flip(tc(2:end)).']; % Let's store this for now.
end

