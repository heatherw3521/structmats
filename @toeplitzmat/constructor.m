function T = constructor(T, varargin)
%CONSTRUCTOR The main TOEPLITZMAT constructor
%   TODO: Detailed explanation goes here

    % CASE: Only first column specified -- assume circulant
    if( numel(varargin) == 1)
        if( isa(varargin{1},'circulantmat') )
            T = toeplitzmat(varargin{1}.tc, varargin{1}.tc);
        else
            T = toeplitzmat(varargin{1},varargin{1});
        end
        return;
    elseif( numel(varargin) > 2)
        error('STRUCTMATS:TOEPLITZMAT:constructor:toomanyargs', ...
                'Too many (>2) arguments passed in to construct Toeplitz matrix');
    end

    % Get the first row and column
    tc = varargin{1};
    tr = varargin{2};
    
    % Error case - bad input pair types
    if( (~isnumeric(tc) && ~islogical(tc)) || ~isa(tr, class(tc)))
        error('STRUCTMATS:TOEPLITZMAT:constructor:inputtypemismatch', ...
                'Inputs pair of types %s and %s is not supported for TOEPLITZMAT.',...
                class(tc), class(tr));
    end
    
    % Flip-flop tc,tr to ensure that they are a column and a row, resp.
    if(size(tc,1) == 1)
        tc = tc.';
    end
    if(size(tr,2) == 1)
        tr = tr.';
    end
    
    % Error case - input size (given col/row are not vectors)
    if( size(tc,2) > 1 && size(tr,1) > 1 )
        error('STRUCTMATS:TOEPLITZMAT:constructor:inputsizemismatch', ...
                'Both input first column and row are not vectors.');
    elseif( size(tc,2) > 1 )
        error('STRUCTMATS:TOEPLITZMAT:constructor:inputsizemismatch', ...
                'Input first column is not a vector.');
    elseif( size(tr,1) > 1 )
        error('STRUCTMATS:TOEPLITZMAT:constructor:inputsizemismatch', ...
                'Input first row is not a vector.');
    end
    
    
    % Error case - inconsistent top-left element
    if(tc(1) ~= tr(1))
        error('STRUCTMATS:TOEPLITZMAT:constructor:firstelementmismatch', ...
                'Given first column and row have different first element.');
    end
    
    % All is good - finish construction
    T.tc = tc;
    T.tr = tr;
end

