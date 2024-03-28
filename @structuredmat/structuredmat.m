classdef (Abstract) structuredmat
    %DIAGONALSTRUCTURE Abstract class for representing structured matrices.
    %   The purpose of this class is to standardize what we require of the 
    %   structured matrices that we implement.
    %   
    %   It is an abstract class, which means that it cannot be
    %   instantiated. However, all of the structured matrix classes we
    %   create (e.g. Toeplitz, Circulant, Hankel, Cauchy, etc.) should
    %   extend this class -- for organizational and practical purposes.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Abstract)

        % full should cast to a regular MATLAB matrix
        r = full(obj1);
        
        % overload the logic and comparison operators
        r = and(obj1,obj2);
        r = or(obj1,obj2);
        r = le(obj1,obj2);
        r = lt(obj1,obj2);
        r = ge(obj1,obj2);
        r = gt(obj1,obj2);
        r = eq(obj1,obj2);
        r = ne(obj1,obj2);
        r = not(obj1);

        % overload the pointwise arithmetic operators
        r = plus(obj1,obj2);
        r = minus(obj1,obj2);
        r = ldivide(obj1,obj2);
        r = rdivide(obj1,obj2);
        r = power(obj1,obj2);
        r = uminus(obj1);
        r = uplus(obj1);

        % overload the matrix arithmetic operators
        r = mtimes(obj1,obj2);
        r = mpower(obj1,obj2);
        r = mldivide(obj1,obj2);
        r = mrdivide(obj1,obj2);

        % overload the subindexing and concatenation
        r = subsref(obj1,obj2);
        r = horzcat(obj1,obj2);
        r = vertcat(obj1,obj2);
        
        % overload the transpose and flip operators
        r = transpose(obj1);
        r = ctranspose(obj1);
        r = flip(obj1);
        r = fliplr(obj1);
    end
end

