classdef hss
    %HSS
    
    properties
        A11 %upper left diag block
        A22 %lower right diag block

        A21 %lower left off diag block
        A12 %upper right off diag block

        sz % [rows cols] of this block; use size(H) from outside the class
        Ir %row indices in H
        Ic %column indices in H
        blocksize %value below which to stop creating new blocks (not necessarily the exact leaf level blocksize)
        levelcount %number of levels
            
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

        Q % QR factor of offdiagonals
        S % RQ factor of diagonals

        minnormQ % QR for underdetermined system. Only stored at root level
        nullorrange % is the minnormQ constructed to span nullspace(A) or range(A^T)
    end

    properties (Hidden, Transient)
        % Stored factorizations (an hssutil.FactorCache handle, root only).
        % H\b, minnorm and tikhonov fill it on their first call and reuse it
        % afterwards. Every assignment to a property that defines the matrix
        % (set methods below) gives the object a new, empty cache, so stale
        % factors are never used; a copy H2 = H shares the factors until one
        % of the two is modified. Transient: factors are not saved to .mat
        % files. clearfactors(H) frees the memory.
        factorcache
    end

    % Property set methods: keep factorcache consistent with the matrix.
    methods
        function obj = set.isroot(obj, v)
            obj.isroot = v;
            if ~isempty(v) && v
                obj.factorcache = hssutil.FactorCache();
            else
                obj.factorcache = [];
            end
        end
        function obj = set.A11(obj, v)
            obj.A11 = v; obj = invalidatefactors(obj);
        end
        function obj = set.A22(obj, v)
            obj.A22 = v; obj = invalidatefactors(obj);
        end
        function obj = set.A21(obj, v)
            obj.A21 = v; obj = invalidatefactors(obj);
        end
        function obj = set.A12(obj, v)
            obj.A12 = v; obj = invalidatefactors(obj);
        end
        function obj = set.sz(obj, v)
            obj.sz = v; obj = invalidatefactors(obj);
        end
        function obj = set.Ir(obj, v)
            obj.Ir = v; obj = invalidatefactors(obj);
        end
        function obj = set.Ic(obj, v)
            obj.Ic = v; obj = invalidatefactors(obj);
        end
        function obj = set.levelcount(obj, v)
            obj.levelcount = v; obj = invalidatefactors(obj);
        end
        function obj = set.level(obj, v)
            obj.level = v; obj = invalidatefactors(obj);
        end
        function obj = set.rowtreeindex(obj, v)
            obj.rowtreeindex = v; obj = invalidatefactors(obj);
        end
        function obj = set.coltreeindex(obj, v)
            obj.coltreeindex = v; obj = invalidatefactors(obj);
        end
        function obj = set.isleaf(obj, v)
            obj.isleaf = v; obj = invalidatefactors(obj);
        end
        function obj = set.isdiag(obj, v)
            obj.isdiag = v; obj = invalidatefactors(obj);
        end
        function obj = set.Z(obj, v)
            obj.Z = v; obj = invalidatefactors(obj);
        end
        function obj = set.Y(obj, v)
            obj.Y = v; obj = invalidatefactors(obj);
        end
        function obj = set.D(obj, v)
            obj.D = v; obj = invalidatefactors(obj);
        end
        function obj = set.lrcomponent(obj, v)
            obj.lrcomponent = v; obj = invalidatefactors(obj);
        end
    end

    methods (Access = private)
        function obj = invalidatefactors(obj)
            % only the root carries a cache; children skip this cheaply
            if ~isempty(obj.factorcache)
                obj.factorcache = hssutil.FactorCache();
            end
        end
    end

    methods (Static, Hidden)
        function obj = loadobj(obj)
            % factorcache is Transient: give a loaded root an empty cache
            if isobject(obj) && ~isempty(obj.isroot) && obj.isroot
                obj.factorcache = hssutil.FactorCache();
            end
            % objects saved before the size -> sz rename (Oct 2026) load
            % without their dimensions: rebuild them from Ir/Ic (or D)
            if isobject(obj) && isempty(obj.sz)
                if ~isempty(obj.Ir) && ~isempty(obj.Ic)
                    obj.sz = [obj.Ir(2)-obj.Ir(1)+1, obj.Ic(2)-obj.Ic(1)+1];
                elseif ~isempty(obj.D)
                    obj.sz = size(obj.D);
                end
            end
        end
    end
    
     methods (Access = public, Static = false )
        function H = hss(varargin)
            % Empty HSS
            if(nargin == 0)
                H.A11 = [];
                H.A22 = [];
                H.A21 = [];
                H.A12 = [];
                H.sz = [];
                H.Ir = [];
                H.Ic = [];
                H.blocksize = [];
                H.levelcount = [];
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

        % function y = mtimes(H1,H2)
        %     y = hss_matvec(H1,H2);
        % end

        function s = subsref(obj,ind)
            switch ind(1).type
                case '()'
                    % H(I,J) and H(K) follow MATLAB's indexing rules: ':',
                    % logical masks, end, and positive integer indices
                    % (any order, repeats allowed). Anything else errors;
                    % previously out-of-range or fractional indices
                    % returned 0 and logical masks were read as indices
                    % 0/1 (audit B03, B04).
                    hss_assertroot(obj, 'indexing');
                    s = hss_index(obj, ind(1).subs);
                    if numel(ind) > 1
                        s = subsref(s, ind(2:end));
                    end
                otherwise
                    s = builtin('subsref',obj,ind);
            end
        end


    end
end

