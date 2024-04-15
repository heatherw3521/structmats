% TESTING DRIVER FOR POINTWISE OPERATIONS ON A STRUCTURED MATRIX

% Test parametersC
M = 500; N = 250;

%% TEST toeplitzmat
sample = @(m,n) SampleToeplitz(m,n);
RunPointwiseTests(sample, M, N);
RunInequalityTests(sample, M, N);

%% TEST circulantmat
sample = @(m,n) SampleCirculant(m);
RunPointwiseTests(sample,M,M); % circulant matrices are square!
RunInequalityTests(sample, M, M);

%% TEST hankelmat
sample = @(m,n) SampleHankel(m,n);
RunPointwiseTests(sample,M,N);
RunInequalityTests(sample, M,N);