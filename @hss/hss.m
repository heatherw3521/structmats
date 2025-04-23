classdef hss
    %HSS
    
    properties
        A11 %upper left diag block
        A22 %lower right diag block

        A21 %lower left off diag block
        A12 %upper right off diag block

        size % size of matrix
        Ir %row indices in H
        Ic %column indices in H
        blocksize %value below which to stop creating new blocks (not necessarily the exact leaf level blocksize)
            
        level % level in tree with 0 being the root
        rowtreeindex % row index relative to the level i am on
        coltreeindex % column index relative to the level i am on

        isleaf %true or false
        isdiag %true or false
        isroot %true or false
        %granular %true of false

        % these contain the indices for the ID rows and columns
        % TODO: ADJUST TO ANY DECOMP NOT JUST ID
        lowrankrows % which rows to use from the full matrix
        lowrankcols % which cols to use from the full matrix
 
        Z % row ID factor
        Y % column ID factor
        D % dense (empty if not a leaf and diagonal)
        lrcomponent % = A(lowrankrows,lowrankcols) (empty if on the diag)
    end
    
     methods (Access = public, Static = false )
        function H = hss(varargin)
            % Empty HSS
            if(nargin == 0)
                H.A11 = [];
                H.A22 = [];
                H.A21 = [];
                H.A12 = [];
                H.size = [];
                H.Ir = [];
                H.Ic = [];
                H.blocksize = [];
                H.level = [];
                H.rowtreeindex = [];
                H.coltreeindex = [];
                H.isleaf = [];
                H.isdiag = [];
                H.isroot = [];
                H.lowrankrows = [];
                H.lowrankcols = [];
                H.Z = [];
                H.Y = [];
                H.D = [];
                H.lrcomponent = [];
                return;
            end

            % Call the main HSS constructor
            H = hss_constructor(H, varargin{:});
        end

        function y = mtimes(H1,H2)
            y = hss_matvec(H1,H2);
        end

        % function sref = subsref(obj,s)
        %       switch s(1).type
        %           case '()'
        %               sref = extract(obj,s);
        %           otherwise 
        %               sref = builtin('subsref',obj,s);
        %       end % switchA
        %   end % function subsref
    end
end

