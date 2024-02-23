classdef circulantmat < toeplitzmat
    %CIRCULANTMAT class for representing Circulant matrices and computing
    %   with them in an efficient manner.
    %
    %   A circulant matrix is defined as a symmetric matrix which is 
    %   constant along its diagonals. This means that it is a symmetric 
    %   Toeplitz matrix. Here is an example of a circulant matrix:
    %       [ 5  4  3  2]
    %       [ 4  5  4  3]
    %   C = [ 3  4  5  4]
    %       [ 2  3  4  5]
    %
    %   The first column of a Circulant matrix determines all
    %   other entries. We can encode a Toeplitz matrix by storing its first row
    %   and first column -- saving a lot of space.
    %
    %   The special structure of these matrices gives rise to highly
    %   efficient algorithms. TODO - explain more
    %
    %   TODO - describe constructors
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        function C = circulantmat(varargin)
            % The main toeplitzmat constructor!
            
            % Error out if no arguments given
            if ( (nargin == 0) )
                error('STRUCTMATS:CIRCULANTMAT:constructor:none', ...
                'Cannot construct a Ciruclant matrix with no inputs');
            end
            % Call the constructor, all the work is done here:
            C = constructor(C, varargin{:});       
        end
    end
end

