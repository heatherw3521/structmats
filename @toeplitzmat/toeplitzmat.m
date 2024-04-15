classdef toeplitzmat < diagonalstructure
    %TOEPLITZMAT class for representing Toeplitz matrices and computing
    %   with them in an efficient manner.
    %
    %   A Toeplitz matrix is defined as a matrix which is constant along
    %   its diagonals. Here is an example of toeplitz matrix:
    %       [ 5  6  7  8]
    %       [ 4  5  6  7]
    %   T = [ 3  4  5  6]
    %       [ 2  3  4  5]
    %       [ 1  2  3  4]
    %
    %   The first column and first row of a Toeplitz matrix determine all
    %   other entries. We can encode a Toeplitz matrix by storing its first row
    %   and first column -- saving a lot of space.
    %
    %   The special structure of these matrices gives rise to highly
    %   efficient algorithms. TODO - explain more
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        function T = toeplitzmat(varargin)
            % The main toeplitzmat constructor!
            
            % Error out if no arguments given
            if ( (nargin == 0) )
                T.tc = [];
                T.tr = [];
            % Call the constructor, all the work is done here:
            else
                T = constructor(T, varargin{:});
            end       
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (GetAccess = public, SetAccess = protected)
        % 'tc', 'tr' are the first column and first row, respectively
        tc 
        tr
    end
end

