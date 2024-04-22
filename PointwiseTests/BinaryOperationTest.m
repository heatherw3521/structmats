classdef BinaryOperationTest < matlab.unittest.TestCase

    % Private / computed properties
    properties (Access=private)
        sampleScalar = @() rand();
        sampleMatrix;
        sampleStructured;
    end

    % Can be specified by user
    properties (Access=public)
        N = 1000; % Can be specified
        M = 1000; % Can be specified
        sample = @(m,n) SampleToeplitz(m,n); % SHOULD BE SPECIFIED BY USER
        binaryFxn = @(x,y) x + y; % SHOULD BE SPECIFIED BY USER
        % Above are just the default values
    end

    methods(TestClassSetup)
        % Shared setup for the entire test class
    end
    
    methods(TestMethodSetup)
        function setSampleMatrix(testCase)
            testCase.sampleMatrix = @() rand(testCase.M, testCase.N);
        end

        function setSampleStructured(testCase)
            testCase.sampleStructured = @() testCase.sample(testCase.M, testCase.N);
        end
    end
    
    methods(Test)
        % Test methods
        
        % scalar on the left
        function scalarLeftPlus(testCase)
            [y, ytrue] = GeneralTest(testCase.sampleScalar, ...
                                   testCase.sampleStructured, ...
                                   testCase.binaryFxn);
            testCase.verifyEqual(full(y), ytrue); % Check answer
        end
        
        % scalar on the right
        function scalarRightPlus(testCase)
            [y, ytrue] = GeneralTest(testCase.sampleStructured, ...
                                   testCase.sampleScalar, ...
                                   testCase.binaryFxn);
            testCase.verifyEqual(full(y), ytrue); % Check answer
        end

        % matrix on the left
        function matrixLeftPlus(testCase)
            [y, ytrue] = GeneralTest(testCase.sampleMatrix, ...
                                   testCase.sampleStructured, ...
                                   testCase.binaryFxn);
            testCase.verifyEqual(full(y), ytrue); % Check answer
        end
        
        % matrix on the right
        function matrixRightPlus(testCase)
            [y, ytrue] = GeneralTest(testCase.sampleStructured, ...
                                   testCase.sampleMatrix, ...
                                   testCase.binaryFxn);
            testCase.verifyEqual(full(y), ytrue); % Check answer
        end
        
        % two structured matrices
        function structuredPlus(testCase)
            [y, ytrue] = GeneralTest(testCase.sampleStructured, ...
                                   testCase.sampleStructured, ...
                                   testCase.binaryFxn);
            testCase.verifyEqual(full(y), ytrue); % Check answer
        end
    end
    
end