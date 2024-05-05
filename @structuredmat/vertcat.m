function C = vertcat(varargin)
%HORZCAT Vertical Concatenation for TOEPLITZMAT objects
 

warning(['Using Naive implementation of VERTCAT.', ...
           'This should be overloaded to take advantage of the structures!']);
C = [full(varargin{1});vertcat(varargin{2:end})];
end

