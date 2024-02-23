% Levent Batakci
% mtimes is blowing my mind...




failed = 0;
passed = 0;
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

disp("PASSED " + passed + "  FAILED " + failed);

function t = runtest(m,n,k)
    
    tc = rand(m,1);
    T = toeplitzmat(tc,[tc(1) rand(1,n-1)]);
    A = toeplitz(T);
    B = rand(n,k);

    if(norm(T*B-A*B,inf) > 1e-12)
        % disp("FAILED!")
        t=false;
    else
        t=true;
        % disp("PASSED!")
    end

end