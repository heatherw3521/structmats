% TESTING DRIVER FOR POINTWISE OPERATIONS ON A STRUCTURED MATRIX

% Test parameters
M = 5000; N = 1250;

%% TEST toeplitzmat
sample = @(m,n) SampleToeplitz(m,n);
RunPointwiseTests(sample, M, N);


%% TEST circulantmat
sample = @(m,n) SampleCirculant(m);
RunPointwiseTests(sample,M,M); % circulant matrices are square!


%% TEST hankelmat
sample = @(m,n) SampleHankel(m,n);
RunPointwiseTests(sample,M,N);