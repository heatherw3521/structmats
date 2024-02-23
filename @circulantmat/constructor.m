function C = constructor(C, varargin)
%CONSTRUCTOR The main CIRCULANTMAT constructor
%   TODO: Detailed explanation goes here

    % Error cases - Wrong number of inputs
    if( numel(varargin) == 0)
        error('STRUCTMATS:CIRCULANTMAT:constructor:missingargs', ...
                'Too few (0) arguments passed in to construct Circulant matrix');
    elseif( numel(varargin) >= 2)
        error('STRUCTMATS:CIRCULANTMAT:constructor:toomanyargs', ...
                'Too many (>1) arguments passed in to construct Circulant matrix');
    end

    % Get the first column
    tc = varargin{1};
    
    % Flip-flop tc,tr to ensure that they are a column and a row, resp.
    if(size(tc,1) == 1)
        tc = tc.';
    end
    
    % Error case - input 'first column' is not a vector
    if( size(tc,2) > 1)
        error('STRUCTMATS:CIRCULANTMAT:constructor:inputsizemismatch', ...
                'Input first column iis not a vector.');
    end
    
    % All is good - finish construction
    C.tc = tc;
    C.tr = C.tc.'; % Let's store this, but we will not compute with it.
    % C.tr = "tc"; % DON'T STORE THIS, JUST ACCESS TC
end

