classdef (Abstract) diagonalstructure
    %DIAGONALSTRUCTURE Abstract class for representing structured matrices
    %   with nice properties along diagonals or anti-diagonals. The big
    %   idea is that applying certain pointwise operations to the matrix
    %   will produce a matrix with the same diagonal structure.
    %
    %   For instance, if all the data of a matrix is encoded by its first
    %   row and column, then pointwise addition with a scalar or matrix of
    %   the same structure can be done more efficiently by simply
    %   performing the operation on the first row/column! In general, this
    %   is O(n+m) instead of O(nm) -- a HUGE saving if n,m are of
    %   comparable size.
    %
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Abstract)
        result = diag_apply(f, fname, ftest, obj1, obj2)
        result = structure_retained(obj, mi, ni); 
        result = full(obj1, obj2);
        
        result = transpose(obj1, obj2);

        result = dsample(obj1, I, J);
    end
end

