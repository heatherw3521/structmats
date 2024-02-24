% Super basic test to try out mtimes for TOEPLITZMAT
% Will be replaced by a corresponding test in a proper
%   testing suite, eventually.



failed = 0; passed = 0;
for m = 10:10:200
    for k = 1:3:50
        t = circtest(m,k);
        if t
            passed = passed + 1;
            % disp("PASSED " + passed)
        else
            failed = failed + 1;
        end
    end
    disp(m)
end
disp("CIRCULANT TEST: passed " + passed + "  failed " + failed);

failed = 0; passed = 0;
for m = 10:10:200
    for n = 5:2:50
        for k = 1:3:50
            t = runtest(m,n,k);
            if t
                passed = passed + 1;
                % disp("PASSED " + passed)
            else
                failed = failed + 1;
            end
        end
    end
    disp(m)
end
disp("TOEPLITZ TEST: passed" + passed + "  failed " + failed);

function t = runtest(m,n,k)
    
    tc = rand(m,1);
    T = toeplitzmat(tc,[tc(1) rand(1,n-1)]);
    A = toeplitz(T);
    B = rand(n,k);

    Y1 = T*B;
    Y2 = A*B;

    if(norm(Y1-Y2,inf) > 1e-11)
        % disp("FAILED!")
        t=false;
    else
        t=true;
        % disp("PASSED!")
    end

end

function t = circtest(m,k)
    
    tc = rand(m,1);
    T = circulantmat(tc);
    A = full(T);
    B = rand(m,k);

    Y1 = T*B;
    Y2 = A*B;

    if(norm(Y1-Y2,inf) > 1e-11)
        % disp("FAILED!")
        t=false;
    else
        t=true;
        % disp("PASSED!")
    end

end