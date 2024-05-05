function RunNormTests(sample, M, N)
%RUNNORMTESTS Temporary testing of matrix norms

A = sample(M,N);
Afull = full(A);

disp("------- running norm tests -------")
% Test the 2-norm
h = norm(A);
htrue = norm(Afull);
if(abs(h-htrue)/abs(htrue) >= 1e-14)
    disp("[FAIL] 2-norm!")
else
    disp("[SUCCESS] 2-norm!")
end

% Test the 1-norm
h = norm(A,1);
htrue = norm(Afull,1);
if(abs(h-htrue)/abs(htrue) >= 1e-14)
    disp("[FAIL] 1-norm!")
else
    disp("[SUCCESS] 1-norm!")
end

% Test the infinity-norm
h = norm(A,'inf');
htrue = norm(Afull,'inf');
if(abs(h-htrue)/abs(htrue) >= 1e-14)
    disp("[FAIL] infinity-norm!")
else
    disp("[SUCCESS] infinity-norm!")
end

% Test the frobenius-norm
h = norm(A,'fro');
htrue = norm(Afull,'fro');
if(abs(h-htrue)/abs(htrue) >= 1e-14)
    disp("[FAIL] frobenius-norm!")
else
    disp("[SUCCESS] frobenius-norm!")
end

disp("------------------")

end

