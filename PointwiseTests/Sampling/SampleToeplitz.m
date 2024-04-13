function T = SampleToeplitz(M,N)
%SAMPLETOEPLITZ Gets a random M-by-N TOEPLITZMAT
    c = rand(M,1);
    r = [c(1) rand(1,N-1)];
    T = toeplitzmat(c,r);
end

