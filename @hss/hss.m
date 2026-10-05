classdef hss
    %HSS  Hierarchically semiseparable (HSS) matrix.
    %   H = hss(A) compresses the dense matrix A (or a function handle
    %   A(I,J) returning entries, with 'sizeA', [m n]) into HSS form by
    %   nested interpolative decompositions. Options: 'blocksize' (leaf
    %   size, default 200), 'tol' (relative ID tolerance, default 1e-12),
    %   'k' (fixed rank instead of a tolerance), 'cutrule', 'decomp'
    %   ('ID' default, 'oldID'; 'SVD' planned).
    %
    %   Supported: H*X, H1*H2, H\b (square or minimum-norm wide solve,
    %   factors kept with H), minnorm, tikhonov, H', H.', H(I,J), full,
    %   size, end, spy, clearfactors.
    %
    %   The matrix is stored as a binary tree: every node holds its four
    %   blocks A11, A12, A21, A22; diagonal leaves hold dense blocks D, and
    %   off-diagonal nodes hold Z*lrcomponent*Y (bases at the leaves,
    %   translation matrices above them).
    
    properties
        A11 % upper left diagonal block (hss node)
        A22 % lower right diagonal block (hss node)

        A21 % lower left off-diagonal block (hss node)
        A12 % upper right off-diagonal block (hss node)

        sz % [rows cols] of this block; use size(H) from outside the class
        Ir % [first last] row of this block in the full matrix
        Ic % [first last] column of this block in the full matrix
        blocksize % leaf size requested at construction (root only)
        levelcount % depth of the tree (level of the leaves)
            
        level % level of this node, 0 at the root
        rowtreeindex % index of this node's row cluster within its level
        coltreeindex % index of this node's column cluster within its level

        isleaf % true for a leaf node
        isdiag % true for a diagonal node
        isroot % true for the root (the whole matrix)

        lowrankrows % skeleton rows selected by the interpolative decomposition
        lowrankcols % skeleton columns selected by the interpolative decomposition
 
        Z % row basis (leaf) or row translation matrix (above the leaves)
        Y % column basis (leaf) or column translation matrix (above the leaves)
        D % dense block of a diagonal leaf (empty otherwise)
        lrcomponent % coupling matrix B of an off-diagonal node: block = Z*B*Y

        Q % unused
        S % unused

        minnormQ % unused
        nullorrange % unused
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
            % objects saved by older versions of the class may lack sz:
            % rebuild it from Ir/Ic (or D)
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
            % hss() is an empty object
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

            H = hss_constructor(H, varargin{:});
        end

        function s = subsref(obj,ind)
            switch ind(1).type
                case '()'
                    % H(I,J) and H(K) follow MATLAB's indexing rules: ':',
                    % logical masks, end, and positive integer indices
                    % (any order, repeats allowed); anything else errors.
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

