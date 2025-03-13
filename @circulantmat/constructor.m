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
tc = reshape(tc,[],1);

% All is good - finish construction
C.tc = tc;
C.tr = [tc(1) flip(tc(2:end)).']; % Let's store this for now.

% Get Sylvester matrix equation quantities
m = length(tc); n = m;
tr = C.tr;

% Form A and B, the circumshift matrices of the right size
% Here, we make the choice of putting '1' in the top right entry
A = speye(m);
A = [A(end,:);A(1:end-1,:)];
C.A = A;
C.B = A;

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
C.U = [u1 u2];
C.V = [v1 v2];
end

