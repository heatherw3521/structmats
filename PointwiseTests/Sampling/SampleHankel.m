function H = SampleHankel(M,N)
%SAMPLEHANKEL Gets a random M-by-N HANKELMAT
    c = rand(M,1);
    r = [c(end) rand(1,N-1)];
    H = hankelmat(c,r);
end

