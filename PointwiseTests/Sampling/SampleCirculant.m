function C = SampleCirculant(M)
%SAMPLECIRCULANT Gets a random M-by-M CIRCULANTMAT
    c = rand(M,1);
    C = circulantmat(c);
end

