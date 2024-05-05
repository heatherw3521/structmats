% TESTING DRIVER FOR POINTWISE OPERATIONS ON A STRUCTURED MATRIX

% Test parametersC
M = 500; N = 500;

%% TEST toeplitzmat
sample = @(m,n) SampleToeplitz(m,n);
RunPointwiseTests(sample, M, N);
RunInequalityTests(sample, M, N);
RunMatrixTests(sample, M, N);
RunNormTests(sample,M,N);
disp("Done testing Toeplitzmat!");

%% TEST circulantmat
sample = @(m,n) SampleCirculant(m);
RunPointwiseTests(sample,M,M); % circulant matrices are square!
RunInequalityTests(sample, M, M);
RunMatrixTests(sample, M, M);
RunNormTests(sample,M,N);
disp("Done testing Circulantmat!");

%% TEST hankelmat
sample = @(m,n) SampleHankel(m,n);
RunPointwiseTests(sample,M,N);
RunInequalityTests(sample, M,N);
RunMatrixTests(sample, M, N);
RunNormTests(sample,M,N);
disp("Done testing Hankelmat!");