function T = constructor(T, varargin)
%CONSTRUCTOR The main TOEPLITZMAT constructor
%   TODO: Detailed explanation goes here

    % Error cases - Wrong number of inputs
    if( numel(varargin) < 2)
        error('STRUCTMATS:TOEPLITZMAT:constructor:missingargs', ...
                'Too few (<2) arguments passed in to construct Toeplitz matrix');
    elseif( numel(varargin) > 2)
        error('STRUCTMATS:TOEPLITZMAT:constructor:toomanyargs', ...
                'Too many (>2) arguments passed in to construct Toeplitz matrix');
    end

    % Get the first row and column
    tc = varargin{1};
    tr = varargin{2};
    
    % Error case - inconsistent top-left element
    if(tc(1) ~= tr(1))
        error('STRUCTMATS:TOEPLITZMAT:constructor:firstelementmismatch', ...
                'Given first column and row have different first element.');
    end

    % All is good - proceed with construction
    T.tc = tc;
    T.tr = tr;
end

