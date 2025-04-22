function H = hss_constructor(H,A,options)

    arguments
        % standard agruments in every case
        H;
        A;
        options.blocksize = 200;
        options.cutrule = @(k) ceil(k/2);
        options.threshold = 1e-12
    
        % if instead of thresholding we wish to prescribe a certain k value
        options.k = 0
    
        % generally we will use interpolative decompositon but this can also be
        % changed
        % alternatives: 'SVD', TODO
        options.decomp = 'ID'
    
    end
    
    blocksize = options.blocksize;
    cutrule = options.cutrule;
    
    % determine how to approximate off-diagonal blocks
    if options.k
        cutoffval = options.k;
        cutofftype = 'k';
    else
        cutoffval = options.threshold;
        cutofftype = 'threshold';
    end
    
    if options.decomp == 'ID'
        decomp = @inter_decomp;
    elseif options.decomp == 'SVD'
        decomp = @svd_decomp;
    else
        error("The given decomposition method must be 'ID or 'SVD'")
    end
    
    
    % toplevel set up
    H.size = size(A);
    H.isroot = true;
    H.isleaf = false;
    H.isdiag = true;
    H.level = 0;
    H.rowtreeindex = 1;
    H.coltreeindex = 1;
    H.Ir = [1,H.size(1)];
    H.Ic = [1,H.size(2)];
    H.blocksize = blocksize; % this keeps track of blocksize at the top level for ease of access, do i need? TODO
    
    % determine where the cuts take place for level 2
    rowcut = cutrule(H.size(1));
    colcut = cutrule(H.size(2));

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
        rowfactorDict,rowindexDict,colfactorDict,colindexDict,blocksize,cutofftype,cutoffval,rr,rc,cutrule,decomp);
    [H.A22,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = recursivestep(A, H.level+1, rr, rc, ...
        2*H.rowtreeindex, 2*H.coltreeindex, ...
        rowfactorDict,rowindexDict,colfactorDict,colindexDict,blocksize,cutofftype,cutoffval,lr,lc,cutrule,decomp);
    

    % off-diag construction at the first level (the remaining are done
    % during the diag recursions)
    [H.A12,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = offdiagconstructor(A,H.level+1,lr, rc, ...
        2*H.rowtreeindex-1,2*H.coltreeindex, ...
        rowfactorDict,rowindexDict,colfactorDict,colindexDict,blocksize,cutofftype,cutoffval,rr,lc,decomp);
    H.A21 = offdiagconstructor(A,H.level+1,rr, lc, ...
        2*H.rowtreeindex,2*H.coltreeindex-1, ...
        rowfactorDict,rowindexDict,colfactorDict,colindexDict,blocksize,cutofftype,cutoffval,lr,rc,decomp);


end

function [H,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = recursivestep(A,level,rows,cols,treerowindex,treecolindex,rowfactorDict,rowindexDict,colfactorDict,colindexDict,blocksize,cutofftype,cutoffval,otherrows,othercols,cutrule,decomp)
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
    H.size = [H.Ir(2)-H.Ir(1)+1,H.Ic(2)-H.Ic(1)+1];
    H.isdiag = true;

    % if size is less than blocksize we are at a leaf node
    if min(H.size)<=blocksize
        H.isleaf = true;
    else
        H.isleaf = false;
        % determine the next split
        rowcut = cutrule(H.size(1));
        colcut = cutrule(H.size(2));
    end

    % if we are at the leaf level
    if H.isleaf
        % store the dense diag block
        H.D = A(H.Ir(1):H.Ir(2),H.Ic(1):H.Ic(2));
   
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
            rowfactorDict,rowindexDict,colfactorDict,colindexDict,blocksize,cutofftype,cutoffval,rr,rc,cutrule,decomp);
        [H.A22,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = recursivestep(A, H.level+1, rr, rc, ...
            2*H.rowtreeindex, 2*H.coltreeindex, ...
            rowfactorDict,rowindexDict,colfactorDict,colindexDict,blocksize,cutofftype,cutoffval,lr,lc,cutrule,decomp);
    
        % construct the offdiags
        [H.A12,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = offdiagconstructor(A,H.level+1,lr, rc, ...
            2*H.rowtreeindex-1,2*H.coltreeindex, ...
            rowfactorDict,rowindexDict,colfactorDict,colindexDict,blocksize,cutofftype,cutoffval,rr,lc,decomp);
        [H.A21,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = offdiagconstructor(A,H.level+1,rr, lc, ...
            2*H.rowtreeindex,2*H.coltreeindex-1, ...
            rowfactorDict,rowindexDict,colfactorDict,colindexDict,blocksize,cutofftype,cutoffval,lr,rc,decomp);
        
    end

end


function [H,rowfactorDict,rowindexDict,colfactorDict,colindexDict] = offdiagconstructor(A,level,rows,cols,treerowindex,treecolindex,rowfactorDict,rowindexDict,colfactorDict,colindexDict,blocksize,cutofftype,cutoffval,otherrows,othercols,decomp)
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
    H.size = [H.Ir(2)-H.Ir(1)+1,H.Ic(2)-H.Ic(1)+1];
    H.isdiag = false;
    if min(H.size)<=blocksize
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
            [z,rows] = decomp(A(H.Ir(1):H.Ir(2),indexsubtraction([1,size(A,2)],othercols)),ctype = cutofftype,cval = cutoffval, orientation = 'rows');
            
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
            [y,cols] = decomp(A(indexsubtraction([1,size(A,1)],otherrows),H.Ic(1):H.Ic(2)),ctype = cutofftype,cval = cutoffval, orientation = 'columns');
    
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
        rows1 = rowindexDict({[H.level+1,2*H.rowtreeindex-1]});
        rows2 = rowindexDict({[H.level+1,2*H.rowtreeindex]});
        
        % decomp on just the selected rows of A
        % TODO: this is specific to ID again :(
        [z,rows] = decomp(A([rows1{1},rows2{1}],indexsubtraction([1,size(A,2)],othercols)),ctype = cutofftype,cval = cutoffval, orientation = 'rows');
        
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
        [y,cols] = decomp(A(indexsubtraction([1,size(A,1)],otherrows),[cols1{1},cols2{1}]),ctype = cutofftype,cval = cutoffval, orientation = 'columns');

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

