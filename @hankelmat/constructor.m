function H = constructor(H, varargin)
%CONSTRUCTOR The main HANKELMAT constructor
%   TODO: Detailed explanation goes here

    % CASE: Only first column specified -- assume circulant
    if( numel(varargin) == 1)
        error('STRUCTMATS:HANKELMAT:constructor:toomanyargs', ...
                'Too few (1) arguments passed in to construct Hankel matrix');
    elseif( numel(varargin) > 2)
        error('STRUCTMATS:HANKELMAT:constructor:toomanyargs', ...
                'Too many (>2) arguments passed in to construct Hankel matrix');
    end

    % Get the first row and column
    hc = varargin{1};
    hr = varargin{2};
    
    % Error case - bad input pair types
    if( (~isnumeric(hc) && ~islogical(hc)) || ~isa(hr, class(hc)))
        error('STRUCTMATS:HANKELMAT:constructor:inputtypemismatch', ...
                'Inputs pair of types %s and %s is not supported for HANKELMAT.',...
                class(hc), class(hr));
    end
    
    % Flip-flop hc,tr to ensure that they are a column and a row, resp.
    if(size(hc,1) == 1)
        hc = hc.';
    end
    if(size(hr,2) == 1)
        hr = hr.';
    end
    
    % Error case - input size (given col/row are not vectors)
    if( size(hc,2) > 1 && size(hr,1) > 1 )
        error('STRUCTMATS:HANKELMAT:constructor:inputsizemismatch', ...
                'Both input first column and row are not vectors.');
    elseif( size(hc,2) > 1 )
        error('STRUCTMATS:HANKELMAT:constructor:inputsizemismatch', ...
                'Input first column is not a vector.');
    elseif( size(hr,1) > 1 )
        error('STRUCTMATS:HANKELMAT:constructor:inputsizemismatch', ...
                'Input first row is not a vector.');
    end
    
    
    % Error case - inconsistent top-left element
    if(hc(end) ~= hr(1))
        error('STRUCTMATS:HANKELMAT:constructor:firstelementmismatch', ...
                'Given first column and row have different bottom left element.');
    end
    
    % All is good - finish construction
    H.hc = hc;
    H.hr = hr;
end

