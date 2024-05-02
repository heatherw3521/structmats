function RunMatrixTests(sample, M, N)
%RUNMATRIXTESTS Summary of this function goes here
%   Detailed explanation goes here

% Set up test basics
tst = BinaryOperationTest;
tst.sample = sample;
tst.M = M; tst.N = N;

% *
tst.binaryFxn = @(x,y) x*y;
run(tst);

% ^ (matrix exponential)
% tst.binaryFxn = @(x,y) x^y;
% run(tst);


%% TODO: Add tests for matrix solves after these are implemented

end

