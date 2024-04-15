function C = horzcat(varargin)
%HORZCAT Horizontal Concatenation for STRUCTUREDMAT objects
  
warning(['Using Naive implementation of HORZCAT.', ...
        'This should be overloaded to take advantage of the structures!']);
C = [full horzcat(varargin{2:end})];
end

