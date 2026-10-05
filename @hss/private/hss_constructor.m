function H = hss_constructor(H,A,options)
%HSS_CONSTRUCTOR  Build an HSS matrix from a dense matrix or an entry function
%   (called by hss). Rows and columns are split together down to a uniform
%   leaf level; every off-diagonal block row and block column is compressed
%   by a nested interpolative decomposition (skeleton rows/columns of the
%   children are reused at the parent).

    arguments
        H;
        A;                 % dense matrix, or function handle A(I,J) with sizeA
        options.blocksize (1,1) double {mustBePositive, mustBeInteger} = 200;
        options.sizeA = size(A);
        options.cutrule = @(k) ceil(k/2);
        options.tol = 1e-12            % relative ID tolerance

        % fixed rank for every off-diagonal block instead of a tolerance (0 = off)
        options.k = 0

        % compression of the off-diagonal blocks: 'ID' (randomized,
        % hssutil.inter_decompv3), 'oldID' (hssutil.inter_decompv2), or
        % 'SVD' (planned, not implemented yet)
        options.decomp = 'ID'

        % passed to hssutil.inter_decompv3; see its header
        options.powerits = 0
        options.escalatemargin = 10
        options.escalatepowerits = 1

    end

    blocksize = options.blocksize;
    % ---- input checks
    % sizeA is only needed when A is a function handle A(I,J) returning
    % entries; for a numeric A it must agree with size(A).
    if isnumeric(A) || islogical(A)
        if ~isequal(double(options.sizeA(:)'), size(A))
            error('hss:sizeA', ['sizeA = %s does not match size(A) = %s; omit sizeA ' ...
                'for a numeric A (it is for function-handle input).'], ...
                mat2str(options.sizeA), mat2str(size(A)));
        end
        if ~all(isfinite(A(:)))
            error('hss:nonFinite', 'A contains NaN or Inf entries.');
        end
    else
        if isa(A, 'function_handle') && isequal(options.sizeA, [1 1])
            error('hss:sizeA', 'For function-handle input pass sizeA = [rows cols].');
        end
        if numel(options.sizeA) ~= 2 || any(options.sizeA < 1) || any(options.sizeA ~= fix(options.sizeA))
            error('hss:sizeA', 'sizeA must be two positive integers [rows cols].');
        end
    end
    options.sizeA = double(options.sizeA(:)');

    % ---- depth of the tree
    % Rows and columns are split together and every leaf sits at the same
    % depth (the solver, discard() and levelup() rely on it). leafLevel is
    % the largest depth at which every leaf, in the worst case (always the
    % floor half), still has both dimensions >= blocksize. A custom cutrule
    % does not change the depth; each of its cuts must lie in 1..k-1
    % (checked by checkcut).
    m = options.sizeA(1);
    n = options.sizeA(2);
    leafLevel = 0;
    while true
        m2 = floor(m/2);
        n2 = floor(n/2);
        if m2 < blocksize || n2 < blocksize
            break
        end
        m = m2;
        n = n2;
        leafLevel = leafLevel + 1;
    end

    if leafLevel == 0
        H.A11 = [];
        H.A22 = [];
        H.A21 = [];
        H.A12 = [];
        H.sz = options.sizeA;
        H.Ir = [];
        H.Ic = [];
        H.blocksize = [];
        H.levelcount = 0;
        H.level = 0;
        H.rowtreeindex = [];
        H.coltreeindex = [];
        H.isleaf = 1;
        H.isdiag = 1;
        H.isroot = 1;
        H.lowrankrows = [];
        H.lowrankcols = [];
        H.Z = [];
        H.Y = [];
        H.D = A(1:options.sizeA(1), 1:options.sizeA(2));
        H.lrcomponent = [];
        return
    end

    cutrule = options.cutrule;

    % determine how to approximate off-diagonal blocks
    if options.k
        cutoffval = options.k;
        cutofftype = 'k';
    else
        cutoffval = options.tol;
        cutofftype = 'threshold';
    end

    if strcmp(options.decomp, 'ID')
        % the sketching options are bound here, so every decomp(...) call
        % below passes only the block, the cutoff and the orientation
        decomp = @(A, varargin) hssutil.inter_decompv3(A, varargin{:}, ...
            powerits = options.powerits, escalatemargin = options.escalatemargin, ...
            escalatepowerits = options.escalatepowerits);
    elseif strcmp(options.decomp, 'oldID')
        decomp = @hssutil.inter_decompv2;
    elseif strcmp(options.decomp, 'SVD')
        error('hss:decomp:notImplemented', "decomp 'SVD' is not implemented yet; use 'ID' or 'oldID'.")
    else
        error('hss:decomp', "decomp must be 'ID', 'oldID' or 'SVD'.")
    end


    % root
    H.sz = options.sizeA;
    H.isroot = true;
    H.isleaf = false;
    H.isdiag = true;
    H.level = 0;
    H.rowtreeindex = 1;
    H.coltreeindex = 1;
    H.Ir = [1,H.sz(1)];
    H.Ic = [1,H.sz(2)];
    H.blocksize = blocksize;
    H.levelcount = leafLevel;

    % split of the root
    rowcut = checkcut(cutrule, H.sz(1));
    colcut = checkcut(cutrule, H.sz(2));

    % interpolation matrices and skeleton indices of every cluster, keyed by
    % [level, index]; the parent's ID reuses its children's skeletons
    rowfactorDict = dictionary();
    rowindexDict = dictionary();
    colfactorDict = dictionary();
    colindexDict = dictionary();


    % row/column ranges of the two children
    lr = [H.Ir(1),H.Ir(1)+rowcut-1];
    rr = [H.Ir(1)+rowcut,H.Ir(2)];
    lc = [H.Ic(1),H.Ic(1)+colcut-1];
    rc = [H.Ic(1)+colcut,H.Ic(2)];

    % diagonal blocks (recursively)
    [H.A11,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = recursivestep(A, H.level+1, lr, lc, ...
        2*H.rowtreeindex-1, 2*H.coltreeindex-1, ...
        rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,rr,rc,cutrule,decomp,options.sizeA);
    [H.A22,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = recursivestep(A, H.level+1, rr, rc, ...
        2*H.rowtreeindex, 2*H.coltreeindex, ...
        rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,lr,lc,cutrule,decomp,options.sizeA);


    % off-diagonal blocks of the root (the deeper ones are built inside the
    % recursion)
    [H.A12,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = offdiagconstructor(A,H.level+1,lr, rc, ...
        2*H.rowtreeindex-1,2*H.coltreeindex, ...
        rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,rr,lc,decomp,options.sizeA);
    H.A21 = offdiagconstructor(A,H.level+1,rr, lc, ...
        2*H.rowtreeindex,2*H.coltreeindex-1, ...
        rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,lr,rc,decomp,options.sizeA);


end

function [H,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = recursivestep(A,level,rows,cols,treerowindex,treecolindex,rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,otherrows,othercols,cutrule,decomp,sizeA)
    % Diagonal node covering rows(1):rows(2), cols(1):cols(2) at the given
    % level, with all its descendants.
    %   A                  dense matrix or entry function
    %   treerowindex/treecolindex   index of the row/column cluster in its level
    %   row/colfactorDict  interpolation matrices per [level, index]
    %   row/colindexDict   skeleton rows/columns per [level, index]
    %   leafLevel          depth of the leaves
    %   cutofftype         'k' (fixed rank) or 'threshold' (relative tolerance)
    %   cutoffval          the rank or the tolerance
    %   otherrows/othercols  row/column range of the sibling block
    %   cutrule            size of the first child of a block of size k
    %   decomp             interpolative decomposition (function handle)
    %   sizeA              size of the full matrix


    H = hss();
    H.isroot = false;
    H.level = level;
    H.Ir = rows;
    H.Ic = cols;
    H.rowtreeindex = treerowindex;
    H.coltreeindex = treecolindex;
    H.sz = [H.Ir(2)-H.Ir(1)+1,H.Ic(2)-H.Ic(1)+1];
    H.isdiag = true;

    if H.level == leafLevel
        H.isleaf = true;
    else
        H.isleaf = false;
        H.levelcount = leafLevel;
        rowcut = checkcut(cutrule, H.sz(1));
        colcut = checkcut(cutrule, H.sz(2));
    end

    if H.isleaf
        % leaf: store the dense diagonal block
        H.D = A(H.Ir(1):H.Ir(2),H.Ic(1):H.Ic(2));
        H.levelcount = leafLevel;
    else
        % row/column ranges of the two children
        lr = [H.Ir(1),H.Ir(1)+rowcut-1];
        rr = [H.Ir(1)+rowcut,H.Ir(2)];
        lc = [H.Ic(1),H.Ic(1)+colcut-1];
        rc = [H.Ic(1)+colcut,H.Ic(2)];

        % diagonal children (recursively)
        [H.A11,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = recursivestep(A, H.level+1, lr, lc, ...
            2*H.rowtreeindex-1, 2*H.coltreeindex-1, ...
            rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,rr,rc,cutrule,decomp,sizeA);
        [H.A22,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = recursivestep(A, H.level+1, rr, rc, ...
            2*H.rowtreeindex, 2*H.coltreeindex, ...
            rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,lr,lc,cutrule,decomp,sizeA);

        % off-diagonal blocks between the two children
        [H.A12,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = offdiagconstructor(A,H.level+1,lr, rc, ...
            2*H.rowtreeindex-1,2*H.coltreeindex, ...
            rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,rr,lc,decomp,sizeA);
        [H.A21,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = offdiagconstructor(A,H.level+1,rr, lc, ...
            2*H.rowtreeindex,2*H.coltreeindex-1, ...
            rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,lr,rc,decomp,sizeA);

    end

end


function [H,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = offdiagconstructor(A,level,rows,cols,treerowindex,treecolindex,rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,otherrows,othercols,decomp,sizeA)
    % Off-diagonal node covering rows(1):rows(2), cols(1):cols(2) (arguments
    % as in recursivestep; otherrows/othercols are the ranges of the diagonal
    % block in the same block row/column, which the IDs exclude). At the
    % leaf level the row basis is an ID of the whole block row of the row
    % cluster (without its diagonal block), the column basis likewise; above
    % the leaves the IDs act on the children's skeleton rows/columns only.
    H = hss();
    H.isroot = false;
    H.level = level;
    H.Ir = rows;
    H.Ic = cols;
    H.rowtreeindex = treerowindex;
    H.coltreeindex = treecolindex;
    H.sz = [H.Ir(2)-H.Ir(1)+1,H.Ic(2)-H.Ic(1)+1];
    H.isdiag = false;
    if H.level == leafLevel
        H.isleaf = true;
    else
        H.isleaf = false;
    end

    if H.isleaf

        % row basis: reuse it if this row cluster was already compressed
        if isConfigured(rowfactorDict) && isKey(rowfactorDict,{[H.level,H.rowtreeindex]})
            z = rowfactorDict({[H.level,H.rowtreeindex]});
            H.Z = z{1};
            rows = rowindexDict({[H.level,H.rowtreeindex]});
            H.lowrankrows = rows{1};
        else
            % row ID of the block row: A(rows,:) = Z * A(skeleton rows,:)
            [z,rows] = decomp(A(H.Ir(1):H.Ir(2),[1:othercols(1)-1,othercols(2)+1:sizeA(2)]),ctype = cutofftype,cval = cutoffval, orientation = 'rows');
            H.Z = z;
            rowfactorDict({[H.level,H.rowtreeindex]}) = {z};

            % skeleton rows as global indices
            H.lowrankrows = rows+H.Ir(1)-1;
            rowindexDict({[H.level,H.rowtreeindex]}) = {rows+H.Ir(1)-1};
        end

        % column basis, likewise
        if isConfigured(colfactorDict) && isKey(colfactorDict,{[H.level,H.coltreeindex]})
            y = colfactorDict({[H.level,H.coltreeindex]});
            H.Y = y{1};
            cols = colindexDict({[H.level,H.coltreeindex]});
            H.lowrankcols = cols{1};
        else
            [y,cols] = decomp(A([1:otherrows(1)-1,otherrows(2)+1:sizeA(1)],H.Ic(1):H.Ic(2)),ctype = cutofftype,cval = cutoffval, orientation = 'columns');
            H.Y = y;
            colfactorDict({[H.level,H.coltreeindex]}) = {y};

            H.lowrankcols = cols+H.Ic(1)-1;
            colindexDict({[H.level,H.coltreeindex]}) = {cols+H.Ic(1)-1};
        end
        % coupling matrix: the skeleton submatrix
        H.lrcomponent = A(H.lowrankrows,H.lowrankcols);

    else
        % above the leaves: IDs of the children's skeleton rows/columns
        rows1 = rowindexDict({[H.level+1,2*H.rowtreeindex-1]});
        rows2 = rowindexDict({[H.level+1,2*H.rowtreeindex]});

        [z,rows] = decomp(A([rows1{1},rows2{1}],[1:othercols(1)-1,othercols(2)+1:sizeA(2)]),ctype = cutofftype,cval = cutoffval, orientation = 'rows');
        H.Z = z;          % row translation matrix
        rowfactorDict({[H.level,H.rowtreeindex]}) = {z};

        oldrows = [rows1{1},rows2{1}];
        H.lowrankrows = oldrows(rows);
        rowindexDict({[H.level,H.rowtreeindex]}) = {H.lowrankrows};

        cols1 = colindexDict({[H.level+1,2*H.coltreeindex-1]});
        cols2 = colindexDict({[H.level+1,2*H.coltreeindex]});

        [y,cols] = decomp(A([1:otherrows(1)-1,otherrows(2)+1:sizeA(1)],[cols1{1},cols2{1}]),ctype = cutofftype,cval = cutoffval, orientation = 'columns');

        H.Y = y;          % column translation matrix
        colfactorDict({[H.level,H.coltreeindex]}) = {y};

        oldcols = [cols1{1},cols2{1}];
        H.lowrankcols = oldcols(cols);
        colindexDict({[H.level,H.coltreeindex]}) = {H.lowrankcols};
        H.lrcomponent = A(H.lowrankrows,H.lowrankcols);
    end
end

function c = checkcut(cutrule, k)
    % cutrule(k) = number of rows (columns) in the first child of a block of
    % size k; must be an integer in 1..k-1
    c = cutrule(k);
    if ~(isscalar(c) && isnumeric(c) && c == fix(c) && c >= 1 && c <= k-1)
        error('hss:cutrule', 'cutrule(%d) returned %s; it must be an integer in 1..%d.', ...
            k, mat2str(c), k-1);
    end
end
