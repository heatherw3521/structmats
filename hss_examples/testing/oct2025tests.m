%% Create A and set parameters
m = 700;
n = 700;
k = 25;
blocksize = 54;

X = linspace(0,1,n);
Y = linspace(0,1,m).' + 0.001;
C = 1 ./ (X-Y);
A = C(1:m,1:n);
function A = genA(m,n)
    X = linspace(0,1,n);
    Y = linspace(0,1,m).' + 0.1;
    C = 1 ./ (X-Y)+eye(m,n);
    A = C(1:m,1:n);
end

%% Construct H
H = hss(A,blocksize = blocksize, k = k);
construct_error = norm(A-full(H))/norm(A);
fprintf('Construction relative error:   %.3e\n', construct_error);

%% Construction + Matvec timing and accuracy tests

% Initialize the number of tests for timing and accuracy
numTests = 5;
cons_timing = zeros(1, numTests);
hssmvec_timing = zeros(1, numTests);
mlmvec_timing = zeros(1, numTests);
mvec_accuracy = zeros(1, numTests);
hssulv_timing = zeros(1, numTests);
mlsolve_timing = zeros(1, numTests);
ulv_accuracy = zeros(1, numTests);

% Measure construction time
for i=1:numTests
    % size
    m = 1e3*i;
    n = 1e3*i;

    % create A, x and b
    A = genA(m,n);
    x_true = randn(n,1);


    % matlab matvec and solver timing
    tic
    b_true = A*x_true;
    mlmvec_timing(i) = toc;
    
    tic
    x_ml = A\b_true;
    mlsolve_timing(i) = toc;



    % specify blocksize and rank for lr approx
    blocksize = m/5;
    k = blocksize/2;

    % time construction
    tic;
    H = hss(A, blocksize = blocksize ,k = k);
    cons_timing(i) = toc;

    % time matvec
    tic
    b = H*x_true;
    hssmvec_timing(i) = toc;
    mvec_accuracy(i) = norm(b - b_true) / norm(b_true);

    % time solver
    tic
    x_ulv = H\b_true;
    hssulv_timing(i) = toc;
    ulv_accuracy(i) = norm(x_ulv-x_true)/norm(x_true);
end

% Plot construction timing
figure
plot(1e3*(1:numTests), cons_timing, '-o');
xlabel('Matrix Size');
ylabel('Construction Time (seconds)');
title('HSS Construction Time');
grid on;

% Plot matvec timing
figure
hold on
plot(1e3*(1:numTests),hssmvec_timing,'-o')
plot(1e3*(1:numTests),mlmvec_timing,'-o')
hold off
xlabel('Matrix Size')
ylabel('Matvec Time (seconds)')
legend('HSS','Matlab')
title('HSS Matvec Timing')

% Plot matvec accuracy
figure
plot(1e3*(1:numTests), mvec_accuracy, '-o');
xlabel('Matrix Size');
ylabel('Matvec Accuracy');
title('HSS Matvec Accuracy');
grid on;


% Plot solver timing
figure
hold on
semilogy(1e3*(1:numTests),hssulv_timing,'-o')
semilogy(1e3*(1:numTests),mlsolve_timing,'-o')
hold off
xlabel('Matrix Size')
ylabel('Solver Time (seconds)')
legend('HSS','Matlab')
title('HSS Solver Timing')

% Plot solver accuracy
figure
semilogy(1e3*(1:numTests), ulv_accuracy, '-o');
xlabel('Matrix Size');
ylabel('Solver Accuracy');
title('HSS Solver Accuracy');
grid on;




%%  Matvec testing

% Random vector and right-hand side
x_true = randn(n,1);
b_true = A * x_true;

% Matvec accuracy
Hx = H * x_true;
matvec_error = norm(Hx - A*x_true) / norm(A*x_true);
fprintf('MatVec relative error:   %.3e\n', matvec_error);

%% 

x_hss = H \ b_true; 
solver_error = norm(x_hss - x_true) / norm(x_true);
residual_error = norm(A*x_hss - b_true) / norm(b_true);

fprintf('Solver relative error:   %.3e\n', solver_error);
fprintf('Residual relative error: %.3e\n', residual_error);


x_ref = A \ b_true;
compare_error = norm(x_hss - x_ref) / norm(x_ref);
fprintf('Compare with dense solve: %.3e\n', compare_error);




%% 
% parameters
m = 400;
n = 400;
k = 7;
blocksize = 43;
% create A
A = genA(m,n);
% create HSS structure for A and visualize
H = hss(A,blocksize = blocksize, k = k); 
%spy(H)
% xtrue and btrue
% xtrue = randn(n,1);
xtrue = ones(n,1);
xtrue(7) = 1;
xtrue(7) = 0;
btrue = A*xtrue;

% matvec error
norm(A*xtrue - H*xtrue)

% HSS ULV solve
x = H\btrue;

% error in solve
norm(x-full(H)\btrue)
norm(x-xtrue)

%% 

N = 5000;
x_all = linspace(-10,10,2*N);
x = x_all(1:2:2*N);
y = x_all(2:2:2*N);
[X,Y] = meshgrid(x,y);
Z = 1 ./ (X-Y) + eye(N);
H = hss(Z);

x = rand(N,1);
norm(H*x-Z*x)
b = Z*x;

norm(x-H\b)