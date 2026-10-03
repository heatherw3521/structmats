function H = hss_constructor(H,A,options)

    arguments
        % standard agruments in every case
        H;
        A;
        % NOTE: was unconstrained (options.blocksize = 200) -- blocksize
        % <= 0 made the leafLevel loop below infinite: it only breaks
        % when m2 or n2 (both >= 0, via floor(.../2)) drops below
        % blocksize, which never happens once blocksize <= 0. Validating
        % here turns that hang into an immediate, clear error instead.
        options.blocksize (1,1) double {mustBePositive, mustBeInteger} = 200;
        options.sizeA = size(A);
        options.cutrule = @(k) ceil(k/2);
        options.tol = 1e-12

        % if instead of thresholding we wish to prescribe a certain k value
        options.k = 0

        % generally we will use interpolative decompositon but this can also be
        % changed
        % alternatives: 'SVD', TODO
        options.decomp = 'ID'

        % passed through to hssutil.inter_decompv3 -- see its own header
        % for what each does. 0 power iterations by default (no change to
        % existing behavior); escalatemargin already defaults sensibly
        % inside inter_decompv3 itself, exposed here for convenience.
        options.powerits = 0
        options.escalatemargin = 10
        options.escalatepowerits = 1

    end

    blocksize = options.blocksize;
    % ---- input checks
    % sizeA is only needed when A is a function handle A(I,J) returning
    % entries. For a numeric A it must agree with size(A): a different
    % value used to build a matrix of the wrong size from part of A
    % (audit B09). NaN/Inf entries used to be accepted and spread through
    % every product (audit B13).
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
    % does not change the depth; each of its cuts is checked to lie in
    % 1..k-1 when it is applied (checkcut), because an unchecked cut such as
    % @(k) k-1 used to create empty blocks and a matrix of the wrong size
    % (audit B08).
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

    % NOTE: was `options.decomp == 'ID'` etc. -- char-array == instead of
    % strcmp, which errors outright ("Arrays have incompatible sizes for
    % this operation") whenever the option string has a different length
    % than whichever literal happened to be compared first, making every
    % branch but the very first unreachable regardless of which one the
    % caller actually asked for.
    if strcmp(options.decomp, 'ID')
        % powerits/escalatemargin baked in via closure rather than threaded
        % through recursivestep/offdiagconstructor's own signatures -- every
        % existing decomp(A_slice, ctype=..., cval=..., orientation=...)
        % call site picks them up automatically through varargin forwarding.
        decomp = @(A, varargin) hssutil.inter_decompv3(A, varargin{:}, ...
            powerits = options.powerits, escalatemargin = options.escalatemargin, ...
            escalatepowerits = options.escalatepowerits);
    elseif strcmp(options.decomp, 'SVD')
        decomp = @svd_decomp;
    elseif strcmp(options.decomp, 'oldID')
        decomp = @hssutil.inter_decompv2;
    else
        error("The given decomposition method must be 'ID', 'oldID', or 'SVD'")
    end


    % toplevel set up
    H.sz = options.sizeA;
    H.isroot = true;
    H.isleaf = false;
    H.isdiag = true;
    H.level = 0;
    H.rowtreeindex = 1;
    H.coltreeindex = 1;
    H.Ir = [1,H.sz(1)];
    H.Ic = [1,H.sz(2)];
    H.blocksize = blocksize; % this keeps track of blocksize at the top level for ease of access, do i need? TODO
    H.levelcount = leafLevel;

    % determine where the cuts take place for level 2
    rowcut = checkcut(cutrule, H.sz(1));
    colcut = checkcut(cutrule, H.sz(2));

    % setup the dictionaries for upward recursion storing decomps at
    % various levels
    rowfactorDict = dictionary();
    rowindexDict = dictionary();
    colfactorDict = dictionary();
    colindexDict = dictionary();


    % indices to recurse into
    lr = [H.Ir(1),H.Ir(1)+rowcut-1];
    rr = [H.Ir(1)+rowcut,H.Ir(2)];
    lc = [H.Ic(1),H.Ic(1)+colcut-1];
    rc = [H.Ic(1)+colcut,H.Ic(2)];

    % diag recursion
    [H.A11,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = recursivestep(A, H.level+1, lr, lc, ...
        2*H.rowtreeindex-1, 2*H.coltreeindex-1, ...
        rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,rr,rc,cutrule,decomp,options.sizeA);
    [H.A22,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = recursivestep(A, H.level+1, rr, rc, ...
        2*H.rowtreeindex, 2*H.coltreeindex, ...
        rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,lr,lc,cutrule,decomp,options.sizeA);


    % off-diag construction at the first level (the remaining are done
    % during the diag recursions)
    [H.A12,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = offdiagconstructor(A,H.level+1,lr, rc, ...
        2*H.rowtreeindex-1,2*H.coltreeindex, ...
        rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,rr,lc,decomp,options.sizeA);
    H.A21 = offdiagconstructor(A,H.level+1,rr, lc, ...
        2*H.rowtreeindex,2*H.coltreeindex-1, ...
        rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,lr,rc,decomp,options.sizeA);


end

function [H,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = recursivestep(A,level,rows,cols,treerowindex,treecolindex,rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,otherrows,othercols,cutrule,decomp,sizeA)
    % inputs:
    % A (mxn array, full matrix)
    % level (integer, current level in recursion with 0 being the root
    % rows (2x1 array, leftmost and rightmost row indices of this block)
    % cols (2x1 array, leftmost and rightmost col indices of this block)
    % treerowindex (integer, row index in the index tree at current level)
    % treecolindex (integer, col index in the index tree at current level)
    % rowfactorDict (dictionary, keys are 2x1 arrays of level and row index, values are the rank k factor for that level and row)
    % colfactorDict (""" but for cols)
    % rowindexDict (dictionary, keys are 2x1 arrays of level and row index, values are k rows computed from an interpolative decomposition of that level and row pairing)
    % colindexDict (""" but for cols)
    % blocksize (integer, size of a leaf block)
    % cutofftype (either 'k' or 'threshold')
    % cutoffval (if cutofftype = k: integer, rank of decomp) (if cutofftype = 'threshold': truncation value at which to stop decomps)
    % otherrows (2x1 array, leftmost and rightmost row indices of the sibling block)
    % othercols (2x1 array, leftmost and rightmost col indices of the sibling block)
    % cutrule (how to split when building the tree)
    % decomp (how to compute the decomposition, either 'ID' or 'svd')


    % create empty hss matrix
    H = hss();

    % properties
    H.isroot = false;
    H.level = level;
    H.Ir = rows;
    H.Ic = cols;
    H.rowtreeindex = treerowindex;
    H.coltreeindex = treecolindex;
    H.sz = [H.Ir(2)-H.Ir(1)+1,H.Ic(2)-H.Ic(1)+1];
    H.isdiag = true;

    % if size is less than blocksize we are at a leaf node
    if H.level == leafLevel
        H.isleaf = true;
    else
        H.isleaf = false;
        H.levelcount = leafLevel;
        % determine the next split
        rowcut = checkcut(cutrule, H.sz(1));
        colcut = checkcut(cutrule, H.sz(2));
    end

    % if we are at the leaf level
    if H.isleaf
        % store the dense diag block
        H.D = A(H.Ir(1):H.Ir(2),H.Ic(1):H.Ic(2));
        H.levelcount = leafLevel;
    % not on the diag? recursion time!!
    else
        % determine the indices for the next level down
        lr = [H.Ir(1),H.Ir(1)+rowcut-1];
        rr = [H.Ir(1)+rowcut,H.Ir(2)];
        lc = [H.Ic(1),H.Ic(1)+colcut-1];
        rc = [H.Ic(1)+colcut,H.Ic(2)];

        % recurse into the diags
        [H.A11,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = recursivestep(A, H.level+1, lr, lc, ...
            2*H.rowtreeindex-1, 2*H.coltreeindex-1, ...
            rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,rr,rc,cutrule,decomp,sizeA);
        [H.A22,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = recursivestep(A, H.level+1, rr, rc, ...
            2*H.rowtreeindex, 2*H.coltreeindex, ...
            rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,lr,lc,cutrule,decomp,sizeA);

        % construct the offdiags
        [H.A12,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = offdiagconstructor(A,H.level+1,lr, rc, ...
            2*H.rowtreeindex-1,2*H.coltreeindex, ...
            rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,rr,lc,decomp,sizeA);
        [H.A21,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = offdiagconstructor(A,H.level+1,rr, lc, ...
            2*H.rowtreeindex,2*H.coltreeindex-1, ...
            rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,lr,rc,decomp,sizeA);

    end

end


function [H,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = offdiagconstructor(A,level,rows,cols,treerowindex,treecolindex,rowfactorDict,rowindexDict,colfactorDict,colindexDict,leafLevel,cutofftype,cutoffval,otherrows,othercols,decomp,sizeA)
    % inputs:
    % A (mxn array, full matrix)
    % level (integer, current level in recursion with 0 being the root
    % rows (2x1 array, leftmost and rightmost row indices of this block)
    % cols (2x1 array, leftmost and rightmost col indices of this block)
    % treerowindex (integer, row index in the index tree at current level)
    % treecolindex (integer, col index in the index tree at current level)
    % rowfactorDict (dictionary, keys are 2x1 arrays of level and row index, values are the rank k factor for that level and row)
    % colfactorDict (""" but for cols)
    % rowindexDict (dictionary, keys are 2x1 arrays of level and row index, values are k rows computed from an interpolative decomposition of that level and row pairing)
    % colindexDict (""" but for cols)
    % blocksize (integer, size of a leaf block)
    % cutofftype (either 'k' or 'threshold')
    % cutoffval (if cutofftype = k: integer, rank of decomp) (if cutofftype = 'threshold': truncation value at which to stop decomps)
    % otherrows (2x1 array, leftmost and rightmost row indices of the sibling block)
    % othercols (2x1 array, leftmost and rightmost col indices of the sibling block)
    % decomp (how to compute the decomposition, either 'ID' or 'svd')


    % creating HSS structure
    H = hss();

    % requisite defs
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

    % if we are at the leaf level
    if H.isleaf

        % if the row decomp has already been computed use it!
        if isConfigured(rowfactorDict) && isKey(rowfactorDict,{[H.level,H.rowtreeindex]})
            z = rowfactorDict({[H.level,H.rowtreeindex]});
            H.Z = z{1};
            rows = rowindexDict({[H.level,H.rowtreeindex]});
            H.lowrankrows = rows{1};
            disp('here rows')
        % otherwise we need to compute it
        else
            % compute decomp
            % TODO: This output is for ID
            % input = z*input(rows,:)
            %[z,rows] = decomp(A(H.Ir(1):H.Ir(2),indexsubtraction([1,size(A,2)],othercols)),ctype = cutofftype,cval = cutoffval, orientation = 'rows');
            [z,rows] = decomp(A(H.Ir(1):H.Ir(2),[1:othercols(1)-1,othercols(2)+1:sizeA(2)]),ctype = cutofftype,cval = cutoffval, orientation = 'rows');
            % store z
            H.Z = z;
            rowfactorDict({[H.level,H.rowtreeindex]}) = {z};

            % store rows
            % adjusting by where in the matrix we are
            H.lowrankrows = rows+H.Ir(1)-1;
            rowindexDict({[H.level,H.rowtreeindex]}) = {rows+H.Ir(1)-1};
        end

        % do the same for columns
        if isConfigured(colfactorDict) && isKey(colfactorDict,{[H.level,H.coltreeindex]})
            y = colfactorDict({[H.level,H.coltreeindex]});
            H.Y = y{1};
            cols = colindexDict({[H.level,H.coltreeindex]});
            H.lowrankcols = cols{1};
            disp('here cols')
        % otherwise we need to compute it
        else
            %[y,cols] = decomp(A(indexsubtraction([1,size(A,1)],otherrows),H.Ic(1):H.Ic(2)),ctype = cutofftype,cval = cutoffval, orientation = 'columns');
            [y,cols] = decomp(A([1:otherrows(1)-1,otherrows(2)+1:sizeA(1)],H.Ic(1):H.Ic(2)),ctype = cutofftype,cval = cutoffval, orientation = 'columns');
            H.Y = y;
            colfactorDict({[H.level,H.coltreeindex]}) = {y};

            % adjusting by where in the matrix we are
            H.lowrankcols = cols+H.Ic(1)-1;
            colindexDict({[H.level,H.coltreeindex]}) = {cols+H.Ic(1)-1};
        end
        % TODO: THIS IS ONLY TRUE FOR ID
        H.lrcomponent = A(H.lowrankrows,H.lowrankcols);

    % not at the leaf level so we have preexisting decomps at a finer level
    else
        % finding the children row decomps
        1;
        rows1 = rowindexDict({[H.level+1,2*H.rowtreeindex-1]});
        rows2 = rowindexDict({[H.level+1,2*H.rowtreeindex]});

        % decomp on just the selected rows of A
        % TODO: this is specific to ID again :(
       % [z,rows] = decomp(A([rows1{1},rows2{1}],indexsubtraction([1,size(A,2)],othercols)),ctype = cutofftype,cval = cutoffval, orientation = 'rows');
        [z,rows] = decomp(A([rows1{1},rows2{1}],[1:othercols(1)-1,othercols(2)+1:sizeA(2)]),ctype = cutofftype,cval = cutoffval, orientation = 'rows');
        % saving z
        H.Z = z;
        rowfactorDict({[H.level,H.rowtreeindex]}) = {z};

        % saving rows
        % adjusting by where in the matrix we are
        oldrows = [rows1{1},rows2{1}];
        H.lowrankrows = oldrows(rows);
        rowindexDict({[H.level,H.rowtreeindex]}) = {H.lowrankrows};

        % columns next
        % decomps exist on finer level
        cols1 = colindexDict({[H.level+1,2*H.coltreeindex-1]});
        cols2 = colindexDict({[H.level+1,2*H.coltreeindex]});

        % TODO: specific to ID
        %[y,cols] = decomp(A(indexsubtraction([1,size(A,1)],otherrows),[cols1{1},cols2{1}]),ctype = cutofftype,cval = cutoffval, orientation = 'columns');
        [y,cols] = decomp(A([1:otherrows(1)-1,otherrows(2)+1:sizeA(1)],[cols1{1},cols2{1}]),ctype = cutofftype,cval = cutoffval, orientation = 'columns');

        H.Y = y;
        colfactorDict({[H.level,H.coltreeindex]}) = {y};

        % adjusting by where in the matrix we are
        oldcols = [cols1{1},cols2{1}];
        H.lowrankcols = oldcols(cols);
        colindexDict({[H.level,H.coltreeindex]}) = {H.lowrankcols};
        H.lrcomponent = A(H.lowrankrows,H.lowrankcols);
    end
end

function Inew = indexsubtraction(I1,I2)
    I1 = I1(1):I1(2);
    I2 = I2(1):I2(2);
    Inew = setdiff(I1,I2);
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
