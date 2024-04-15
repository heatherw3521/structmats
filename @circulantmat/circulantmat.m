classdef circulantmat < toeplitzmat
    %CIRCULANTMAT class for representing Circulant matrices and computing
    %   with them in an efficient manner.
    %
    %   A circulant matrix is a special type of Toeplitz matrix which is 
    %   determined by its first column (or first row). The key property of
    %   a circulant matrix is that column j+1 is a circular shift down of
    %   column j. Here is an example of a circulant matrix:
    %       [ 5  2  3  4]
    %       [ 4  5  2  3]
    %   C = [ 3  4  5  2]
    %       [ 2  3  4  5]
    %
    %   The first column of a Circulant matrix determines all
    %   other entries. We can encode a Circulant matrix by storing its first row
    %   and first column -- saving a lot of space.
    %
    %   The special structure of these matrices gives rise to highly
    %   efficient algorithms. TODO - explain more
    %
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        function C = circulantmat(varargin)
            % The main toeplitzmat constructor!
            
            % Initialize empty TOEPLITZMAT
            if ( (nargin == 0) )
                error('STRUCTMATS:CIRCULANTMAT:constructor:none', ...
                'Cannot construct a Ciruclant matrix with no inputs');
            end
            % Call the constructor, all the work is done here:
            C = constructor(C, varargin{:});       
        end
    end
end

