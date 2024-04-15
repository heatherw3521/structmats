classdef vandermat < structuredmat
    %VANDERMAT class for representing Vandermonde matrices and computing
    %   with them in an efficient manner.
    %
    %   A Vandermonde matrix is defined by its rows forming geometric 
    %   progressions Vij = (V2j)^i. Here is a small example
    %       [ 1  2  4   8  ]
    %       [ 1  3  9   27 ]
    %   V = [ 1  4  16  64 ]
    %       [ 1  5  25  125]
    %
    %   If the Vandermonde matrix's dimension is known, then the second
    %   column determines the rest of the matrix. For efficiency purposes,
    %   this will be how we store the matrix.
    % 
    %   Vandermonde matrices are often used in problems invovling
    %   polynomials. Indeed, given a vector of coefficients c, the
    %   matrix-vector product Vc computes
    %
    %        [ 1  x1  (x1)^2] [c0]   [p(x1)]
    %   Vc = [ 1  x2  (x2)^2] [c1] = [p(x2)] ,
    %        [ 1  x3  (x2)^2] [c2]   [p(x3)]
    %   where p(x) is the polynomial defined by
    %       p(x) = c0 + c1*x + c2*x^2.
    %   
    %   The example above shows us that if V is invertible (which is
    %   guaranteed by the xj's being distinct and V being square), we can
    %   use Vandermonde matrices to solve for the coefficients of a
    %   polynomial. Moreover, we can also use Vandermonde matrices to solve
    %   for polynomial best-fits by instead using least-squares on
    %   overdetermined system. It should be noted that a Vandermonde 
    %   matrix need not be square.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        function V = vandermat(varargin)
            % The main vandermat constructor!
            
            error(['VANDERMAT is a work in progress! All operations are ',...
                    'currently implemented via casting to a MATLAB matrix.'])

            % Error out if no arguments given
            if ( (nargin == 0) )
                error('Cannot construct a Vandermonde matrix with no inputs!')
            % Call the constructor, all the work is done here:
            else
                V = constructor(V, varargin{:});
            end       
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (GetAccess = public, SetAccess = protected)
        % 'x', 'n' are the second column and number of columns, respectively
        x
        n
    end
end

