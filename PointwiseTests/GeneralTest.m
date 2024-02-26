function [y, ytrue] = GeneralTest(sampler1, sampler2, fxn)
%GENERALTEST helper function that has the structure of just about any test
%   we'd like to run. Most operations take we'd like to test are binary,
%   and the two different "sampler"s in the input give us all the
%   flexibility we need.

    A = sampler1();
    B = sampler2();

    y = fxn(A,B);
    ytrue = fxn(full(A),full(B));
end

