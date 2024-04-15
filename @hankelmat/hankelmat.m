classdef hankelmat < diagonalstructure
    %HANKELMAT class for representing Hankel matrices and computing
    %   with them in an efficient manner.
    %
    %   A HANKEL matrix is defined as a matrix which is constant along
    %   its antidiagonals. Here is an example of hankel matrix:
    %       [ 1  2  3  4]
    %       [ 2  3  4  5]
    %   T = [ 3  4  5  6]
    %       [ 4  5  6  7]
    %       [ 5  6  7  8]
    %
    %   The first column and last row of a Hankel matrix determine all
    %   other entries. We can encode a Hankel matrix by storing its last row
    %   and first column -- saving a lot of space.
    %
    %   The special structure of these matrices gives rise to highly
    %   efficient algorithms. TODO - explain more
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        function H = hankelmat(varargin)
            % The main toeplitzmat constructor!
            
            % Error out if no arguments given
            if ( (nargin == 0) )
                error(['STRUCTMATS:HANKELMAT:constructor:toofewargs', newline, ...
                   'Too few (0) arguments passed in to construct Hankel matrix']);  
            % Call the constructor, all the work is done here:
            else
                H = constructor(H, varargin{:});
            end       
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (GetAccess = public, SetAccess = protected)
        % 'tc', 'tr' are the first column and first row, respectively
        hc 
        hr
    end
end

