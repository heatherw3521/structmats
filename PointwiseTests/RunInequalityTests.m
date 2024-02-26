function RunInequalityTests(sample, M, N)
%RUNPOINTWISETESTS Summary of this function goes here
%   Detailed explanation goes here

% Set up test basics
tst = BinaryOperationTest;
tst.sample = sample;
tst.M = M; tst.N = N;

% >
tst.binaryFxn = @(x,y) x>y;
run(tst);

% >=
tst.binaryFxn = @(x,y) x>=y;
run(tst);

% <
tst.binaryFxn = @(x,y) x<y;
run(tst);

% <=
tst.binaryFxn = @(x,y) x<=y;
run(tst);

%% TO DO: EDIT THESE TESTS TO USE A SMARTER RANDOM SAMPLING
% ~=
tst.binaryFxn = @(x,y) x~=y;
run(tst);

% ==
tst.binaryFxn = @(x,y) x~=y;
run(tst);

end

