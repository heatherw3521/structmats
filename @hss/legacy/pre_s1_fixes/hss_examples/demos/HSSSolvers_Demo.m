%% 1. Lets create a couple matrices that accept HSS structure
%% kernel instead of matrix
N = 2000;
A = Hankel_kernel(N);
% fullA = A(1:N,1:N);
% Create the HSS representation for the initial matrix A
tic
HZ_A = hss(A,blocksize = 100, sizeA = [N,N]);
toc

b = rand(N,1);

% tic
% xmat = fullA*b;
% toc

tic
xhss = HZ_A*b;
toc

% norm(xmat-xhss)

%% Crash-proof scaling test for Hankel kernel HSS
clear; clc;

Ns = round(logspace(3,5,8));   % 1e3 ... 1e5
blocksize = 100;

% Hard safety limits (tune once, then forget)
MAX_FULL_N = 15000;      % absolute cutoff (VERY important)
MEMORY_FRAC = 0.30;      % use at most 30% of available memory

results = struct( ...
    'N',        Ns, ...
    'hss_build',nan(size(Ns)), ...
    'hss_mv',   nan(size(Ns)), ...
    'full_mv',  nan(size(Ns)), ...
    'full_ok',  false(size(Ns)) ...
);

for k = 1:length(Ns)
    N = Ns(k);
    fprintf('\n==============================\n');
    fprintf('N = %d\n', N);
    fprintf('==============================\n');

    % Kernel operator
    A = Hankel_kernel(N);

    % Build HSS (kernel)
    fprintf('Building HSS (kernel)...\n');
    tic
    HZ_A = hss(A, blocksize = blocksize, sizeA = [N,N]);
    results.hss_build(k) = toc;

    % HSS matvec
    b = rand(N,1);
    tic
    xhss = HZ_A * b;
    results.hss_mv(k) = toc;
    fprintf('HSS matvec time: %.4f s\n', results.hss_mv(k));

    % Decide *beforehand* whether full matrix is allowed
    do_full = false;

    if N <= MAX_FULL_N
        try
            m = memory;
            bytes_needed = 16 * N^2;   % complex double
            if bytes_needed < MEMORY_FRAC * m.MemAvailableAllArrays
                do_full = true;
            else
                fprintf('⚠️  Full matrix skipped (memory estimate too large).\n');
            end
        catch
            % memory() unavailable → rely on hard cutoff
            do_full = true;
        end
    else
        fprintf('⚠️  Full matrix skipped (N > %d hard limit).\n', MAX_FULL_N);
    end

    % Full matrix matvec (only if explicitly allowed)
    if do_full
        fprintf('Forming full matrix (SAFE)...\n');
        fullA = A(1:N,1:N);   % now this is safe

        tic
        xfull = fullA * b;
        results.full_mv(k) = toc;
        results.full_ok(k) = true;

        fprintf('Full matvec time: %.4f s\n', results.full_mv(k));

        % Optional accuracy check
        relerr = norm(xfull - xhss) / norm(xfull);
        fprintf('Relative error: %.2e\n', relerr);
    end

    clear HZ_A fullA A
end


figure; hold on;
loglog(results.N, results.hss_mv, 'o-', 'LineWidth', 2);

idx = results.full_ok;
loglog(results.N(idx), results.full_mv(idx), 's--', 'LineWidth', 2);

xlabel('N');
ylabel('Matvec time (s)');
legend('HSS (kernel)', 'Full matrix', 'Location','NorthWest');
grid on;

%% Hankel kernel HSS scaling benchmark v2
clear; clc;

% Parameters
Ns = round(logspace(3,5,12));   % 12 logspaced sizes
N_full_limit = 20000;          % do NOT form full matrix above this
nreps = 20;                     % number of matvec repetitions
save_hss_structs = true;      % set true if you want to save HSS objects

% Storage
hss_build_time = zeros(size(Ns));
hss_mv_time    = zeros(size(Ns));
full_mv_time   = nan(size(Ns));

if save_hss_structs
    HSS_cells = cell(size(Ns));
end

% Main loop
for k = 1:length(Ns)

    N = Ns(k);
    fprintf('\n====================================\n');
    fprintf('N = %d\n', N);
    fprintf('====================================\n');

    % Kernel
    A = Hankel_kernel(N);

    % --- HSS Construction ---
    blocksize = min(max(100, round(0.5*sqrt(N))),2000);
    fprintf('Building HSS...\n');
    tic
    HZ_A = hss(A, blocksize = blocksize, sizeA = [N,N]);
    hss_build_time(k) = toc;
    fprintf('HSS build time: %.3f s\n', hss_build_time(k));

    if save_hss_structs
        HSS_cells{k} = HZ_A;
    end

    % --- HSS Matvec (averaged) ---
    b = rand(N,1);

    tic
    for r = 1:nreps
        x = HZ_A * b;
    end
    total1 = toc;

    hss_mv_time(k) = total1 / nreps;
    fprintf('Avg HSS matvec time: %.6f s\n', hss_mv_time(k));

    % --- Full matrix (only small N) ---
    if N <= N_full_limit

        fprintf('Forming full matrix...\n');
        fullA = A(1:N,1:N);

        tic
        for r = 1:nreps
            xfull = fullA * b;
        end
        total1 = toc;

        full_mv_time(k) = total1 / nreps;
        fprintf('Avg full matvec time: %.6f s\n', full_mv_time(k));

        clear fullA
    else
        fprintf('Skipping full matrix (N > %d)\n', N_full_limit);
    end

    clear HZ_A A
end

% =============================
% Save Results
% =============================

results_table = table(Ns(:), ...
                      hss_build_time(:), ...
                      hss_mv_time(:), ...
                      full_mv_time(:), ...
    'VariableNames', {'N','HSS_build','HSS_matvec','Full_matvec'});

% Save as text
writetable(results_table, 'hss_scaling_results.txt', 'Delimiter','\t');

% Save as MAT file (includes everything)
if save_hss_structs
    save('hss_scaling_results_v2.mat', ...
         'results_table', ...
         'HSS_cells');
else
    save('hss_scaling_results_v2.mat', ...
         'results_table');
end

fprintf('\nBenchmark complete. Results saved.\n');

%% plot matvec timings
load hss_scaling_results_v2.mat

N  = results_table.N;
tH = results_table.HSS_matvec;

figure; hold on;

% --- Actual HSS matvec ---
loglog(N, tH, 'o-', 'LineWidth', 2);

% --- Dense (if available) ---
idx = ~isnan(results_table.Full_matvec);
if any(idx)
    loglog(N(idx), results_table.Full_matvec(idx), ...
           's--', 'LineWidth', 2);
end

% --- N log N reference curve ---
ref = N .* log(N);

% Scale reference to match last valid HSS point
c = tH(end) / ref(end);
ref_scaled = c * ref;

loglog(N, ref_scaled, 'k--', 'LineWidth', 2);

grid on;
xlabel('N');
ylabel('Average matvec time (s)');
legend('HSS', 'Full', 'c N log N', 'Location','NorthWest');

title('HSS Matvec Timing');



%% ULV Solve Benchmark (from loaded HSS matrices)
clear; clc;

% Load HSS data
load hss_scaling_results_v2.mat   % loads results_table + HSS_cells
%Ns = results_table.N;
Ns = round(logspace(3,5,12));   % 12 logspaced sizes

N_full_limit = 20000;
nreps = 5;

% Load existing results if available
if isfile('ulv_scaling_results.mat')
    load('ulv_scaling_results.mat','solve_table');

    % Remove Dense_fact column if it exists
    if any(strcmp('Dense_fact', solve_table.Properties.VariableNames))
        solve_table.Dense_fact = [];
    end

else
    solve_table = table([],[],[], ...
        'VariableNames', {'N','ULV_solve','Dense_solve'});
end

% Determine which N need to run
needs_run = false(size(Ns));

for i = 1:length(Ns)
    N = Ns(i);
    row = find(solve_table.N == N);

    if isempty(row)
        needs_run(i) = true;              % never computed
    elseif solve_table.ULV_solve(row) == 0
        needs_run(i) = true;              % incomplete timing
    end
end

Ns_to_run = Ns(needs_run);
fprintf('\nNs remaining to run:\n');
disp(Ns_to_run)

% Run only missing sizes
for k = 1:length(Ns_to_run)

    N = Ns_to_run(k);
    idx = find(Ns == N);

    fprintf('\n====================================\n');
    fprintf('Running N = %d\n', N);
    fprintf('====================================\n');

    HZ_A = HSS_cells{idx};

    % --- ULV Solve (averaged) ---
    total1 = 0;
    for r = 1:nreps
        b = rand(N,1);
        tic
        x = HZ_A \ b;
        total1 = total1 + toc;
    end
    ulv_time = total1 / nreps;
    fprintf('Avg ULV solve time: %.6f s\n', ulv_time);

    % --- Dense solve (small N only) ---
    dense_time = NaN;
    if N <= N_full_limit
        fprintf('Reconstructing dense matrix...\n');
        A_dense = full(HZ_A);

        total1 = 0;
        for r = 1:nreps
            b = rand(N,1);
            tic
            xfull = A_dense \ b;
            total1 = total1 + toc;
        end
        dense_time = total1 / nreps;

        fprintf('Avg dense solve time: %.6f s\n', dense_time);

        clear A_dense
    else
        fprintf('Skipping dense solve (N > %d)\n', N_full_limit);
    end

    % Append new row
    new_row = table(N, ulv_time, dense_time, ...
        'VariableNames', {'N','ULV_solve','Dense_solve'});

    solve_table = [solve_table; new_row];

    % Sort and save progress (crash-safe)
    solve_table = sortrows(solve_table,'N');
    save('ulv_scaling_results.mat','solve_table');
    writetable(solve_table,'hss_ulv_scaling_results.txt','Delimiter','\t');

    clear HZ_A
end

fprintf('\nULV benchmark complete. Results appended and saved.\n');

%% ================= Plot =================
load ulv_scaling_results.mat
N = solve_table.N;

figure;
loglog(N, solve_table.ULV_solve, 'o-','LineWidth',2);
hold on;
idx = ~isnan(solve_table.Dense_solve);
if any(idx)
    loglog(N(idx), solve_table.Dense_solve(idx), 's--','LineWidth',2);
end

% N log N reference
ref = N .* log(N).^2;
c = solve_table.ULV_solve(end-7)/ref(end-7);
loglog(N, c*ref, 'k--','LineWidth',2);

grid on;
xlabel('N');
ylabel('Average solve time (s)');
legend('ULV','Dense','c N log^2 N','Location','NorthWest');
title('ULV Solve Scaling');

%% ULV Solve Benchmark (rerun ULV only, keep dense times)
clear; clc;

% Load HSS data
load hss_scaling_results_v2.mat   % loads results_table + HSS_cells
Ns = results_table.N;
nreps = 5;

% Load existing solve table (must exist)
if ~isfile('ulv_scaling_results.mat')
    error('ulv_scaling_results.mat not found.');
end
load('ulv_scaling_results.mat','solve_table');

% Remove Dense_fact column if it exists
if any(strcmp('Dense_fact', solve_table.Properties.VariableNames))
    solve_table.Dense_fact = [];
end

% Add ULV_factored column if it doesn't exist
if ~any(strcmp('ULV_factored', solve_table.Properties.VariableNames))
    solve_table.ULV_factored = nan(height(solve_table), 1);
end

% Rerun ULV for ALL N
for k = 1:length(Ns)
    N = Ns(k);
    idx = find(Ns == N);
    
    fprintf('\n====================================\n');
    fprintf('Recomputing ULV for N = %d\n', N);
    fprintf('====================================\n');
    
    HZ_A = HSS_cells{idx};

    % --- ULV Solve (averaged) ---
    total1 = 0;
    total2 = 0;
    for r = 1:nreps
        b = rand(N,1);
        tic
        [x,HZ_A_stored] = HZ_A \ b;
        total1 = total1 + toc;
        tic
        x2 = HZ_A_stored\b;
        total2 = total2+toc;
    end
    ulv_time = total1 / nreps;
    ulv_factor_time = total2/nreps;
    
    fprintf('Avg ULV solve time: %.6f s\n', ulv_time);
    fprintf('Avg ULV solve time with stored factors: %.6f s\n', ulv_factor_time);

    % Overwrite ULV_solve and ULV_factored but KEEP Dense_solve
    row = find(solve_table.N == N);
    if isempty(row)
        % If somehow missing, append with existing dense = NaN
        new_row = table(N, ulv_time, ulv_factor_time, NaN, ...
            'VariableNames', {'N','ULV_solve','ULV_factored','Dense_solve'});
        solve_table = [solve_table; new_row];
    else
        solve_table.ULV_solve(row) = ulv_time;
        solve_table.ULV_factored(row) = ulv_factor_time;
    end
    
    solve_table = sortrows(solve_table,'N');

    % Save progress after each N (crash-safe)
    save('ulv_scaling_results.mat','solve_table');
    writetable(solve_table,'hss_ulv_scaling_results.txt','Delimiter','\t');
    
    clear HZ_A HZ_A_stored
end

fprintf('\nULV recomputation complete.\n');

%% ================= Plot =================
load ulv_scaling_resultsv2.mat
N = solve_table.N;

figure;
loglog(N, solve_table.ULV_solve, 'o-','LineWidth',2);
hold on;
idx = ~isnan(solve_table.Dense_solve);
if any(idx)
    loglog(N(idx), solve_table.Dense_solve(idx), 's--','LineWidth',2);
end

% N log N reference
ref = N .* (log(N)).^2;
c = solve_table.ULV_solve(end)/ref(end);
loglog(N, c*ref, 'k--','LineWidth',2);

grid on;
xlabel('N');
ylabel('Average solve time (s)');
legend('ULV','Dense','c N log^2 N','Location','NorthWest');
title('ULV Solve Timing');


%% a) Square 
N = 100;
Z = Hankel_mat(N);
x = rand(N,3);


% create the HSS representation
tic
HZ1 = hss(Z,blocksize = 10)
t = toc;
disp("Size = " + N + " x " + N)
disp("HSS construction: " + t + " seconds")

% test the matvec
tic
Z*x;
t = toc;
disp("Matlab matvec: " + t + " seconds")

tic
HZ1*x;
t = toc;
disp("HSS matvec: " + t + " seconds")

disp("Solution Relative Error: " + norm(Z*x-HZ1*x)/norm(Z*x))
figure
spy(HZ1)

%% lets up the size of the square
N = 10000;
Z = Hankel_mat(N);
x = rand(N,1);
disp("----------")
disp("Size = " + N + " x " + N)
% create the HSS representation
tic
HZ = hss(Z);
t = toc;
disp("HSS construction: " + t + " seconds")

% test the matvec
tic
Z*x;
t = toc;
disp("Matlab matvec: " + t + " seconds")

tic
HZ*x;
t = toc;
disp("HSS matvec: " + t + " seconds")

disp("Solution Relative Error: " + norm(Z*x-HZ*x)/norm(Z*x))

%% We can be a bit more precise in construction
% create the HSS representation
tic
HZ2 = hss(Z,blocksize = 100);
t = toc;
disp("----------")
disp("HSS construction: " + t + " seconds")

% test the matvec
tic
Z*x;
t = toc;
disp("Matlab matvec: " + t + " seconds")

tic
HZ2*x;
t = toc;
disp("HSS matvec: " + t + " seconds")

disp("Solution Relative Error: " + norm(Z*x-HZ2*x)/norm(Z*x))
% figure
% spy(HZ)
% figure
% spy(HZ2)
%% Next, the ULV solver!
N = 5000;
Z5k = Hankel_mat(N);
b = rand(N,1);
disp("----------")
disp("Size = " + N + " x " + N)
tic
HZ5k = hss(Z5k,blocksize = 500);
t = toc;
disp("HSS construction: " + t + " seconds")

tic
xtrue = Z5k\b;
t = toc;
disp("Matlab backslash: " + t + " seconds")


tic
xhss = HZ5k\b;
t = toc;
disp("HSS backslash: " + t + " seconds")

disp("Solution Relative Error: " + norm(xhss-xtrue)/norm(xtrue))
figure
spy(HZ5k)

%% How much does storing factors help us?
count = 50;
bs = rand(N,count);
mlt = 0;
hsst_save = 0;
hsst_nosave = 0;
for i=1:count
    b = bs(:,i);
    tic
    xtrue = Z5k\b;
    mlt = toc/count + mlt;
    tic
    xhss = HZ5k\b;
    hsst_save = toc/count + hsst_save;

    Htemp = hss(Z5k,blocksize = 500,tol = 1e-10);
    tic
    xhss = Htemp\b;
    hsst_nosave = toc/count + hsst_nosave;
end



%% And what if we take it back to N=10000?
N=10000;
b = rand(N,1);
disp("----------")
disp("Size = " + N + " x " + N)

tic
xtrue = Z\b;
t = toc;
disp("Matlab backslash: " + t + " seconds")

tic
xhss = HZ\b;
disp("HSS backslash: " + toc + " seconds")

disp("Solution Relative Error: " + norm(xhss-xtrue)/norm(xtrue))


%% Squares seem to be doing ok. I'm a fan of rectangles though, what gives?

M = 2327;
N = 5392;
k = 50;
blocksize = 200;
[Z,HZ_rect,time] = rectsampler(M, N, k, blocksize);
disp("----------")
disp("Size = " + M + " x " + N)
disp("HSS construction: " + time + " seconds")
figure
spy(HZ_rect)

%%
H = hss(randn(M,N),blocksize = 200,tol = 0.8);
spy(H)

%% 
x = rand(N,1);
b = Z*x;


tic
xmat = Z\b;
toc
norm(Z*xmat-b)
norm(xmat)
disp('----')

tic
xulv = HZ_rect\b;
toc
norm(Z*xulv-b)
norm(xulv)
disp('----')
% 
% tic
% xpcg = ud_normeqs_pcg(HZ_rect,b);
% toc
% norm(Z*xpcg-b)
% norm(xpcg)
% disp('----')

tic
xtr = ud_tr_projgrad(HZ_rect,b,xulv);
toc
norm(Z*xtr-b)
norm(xtr)
disp('----')



%% matvec rectangle
x = randn(N,100);

tic
Z*x;
t = toc;
disp("Matlab matvec: " + t + " seconds")

tic
HZ_rect*x;
t = toc;
disp("HSS matvec: " + t + " seconds")


disp("Solution Relative Error: " + norm(Z*x-HZ_rect*x)/norm(Z*x))
%% Lets solve
x = rand(N,1);
b = Z*x;

tic
xtrue = lsqminnorm(Z,b);
t = toc;
disp("Matlab lsminnorm: " + t + " seconds")

tic
[xhss,HZ_rect] = HZ_rect\b; %#ok<RHSFN>
t = toc;
disp("HSS lsminnorm: " + t + " seconds")


tic
fun = @(x, mode) lsqr_wrapper(x, mode, HZ_rect);
[x_hss_lsqr, flag] = lsqr(fun, b, 0, 1000);
t = toc;
disp("HSS lsqrminnorm: " + t + " seconds")


disp("Solution Relative Error: " + norm(xhss-xtrue)/norm(xtrue))

norm(Z*xhss-b)
norm(Z*xtrue-b)
norm(Z*x_hss_lsqr-b)
norm(xhss)
norm(xtrue)
norm(x_hss_lsqr)

%%

HZt_rect = HZ_rect.';

HHt = full(HZ_rect)*full(HZt_rect);
% HHt = hss(ZZt);

lambda = HHt\(-b);

x = -(HZt_rect*lambda);
norm(x)
norm(Z*x-b)

lambda = (Z*Z')\(-b);
x = -Z'*lambda;
norm(x)
norm(Z*x-b)

%% Now the factors are saved does that help?
disp("----------")
disp("Following one solve, we store the QR factor")

x = rand(N,1);
b = Z*x;

tic
xtrue = lsqminnorm(Z,b);
t = toc;
disp("Matlab lsminnorm: " + t + " seconds")

tic
xhss = HZ_rect\b; %#ok<RHSFN>
t = toc;
disp("HSS lsminnorm: " + t + " seconds")

disp("Solution Relative Error: " + norm(xhss-xtrue)/norm(xtrue))


%% obsolete

% x0 = xhss;
% H = HZ_rect;
% iters = 500;
% 
% timing = zeros(iters,1);
% error = zeros(iters,1);
% xnorm = zeros(iters,1);
% 
% % for i=1:iters
% %     tic
% %     xmin1 = minnormv1(HZ_rect,xhss,b,i);
% %     timing(i) = toc;
% %     error(i) = norm(Z*xmin1-b);
% %     xnorm(i) = norm(xmin1);
% % end
% N = size(x0,1);
% Null = zeros(N,iters);
% for i=1:iters
%     v = randn(N,1);
%     bk = b+H*v; %nlogn
%     xk = H\bk; %
%     z = xk-x0-v;
%     Null(:,i) = z;
% end
% 
% for i=1:iters
%     [Q,~] = qr(Null(:,1:i),'econ');
%     xmin = x0-Q*Q'*x0;
%     if norm(xmin-xtrue)<1e-4
%         closeidx = i;
%     end
%     xnorm(i) = norm(xmin);
% end


%% obsolete
% plot(1:iters,timing)

% plot(1:iters,xnorm,'-',LineWidth=1)
% hold on
% plot(1:iters,repelem(norm(xtrue),iters),'--',Color='r',LineWidth=0.1)
% xline(closeidx)
% xlabel('Iterations')
% ylabel('norm')
% legend('norm(x)','norm(xtrue)')

% tic
% xmin2 = minnormv2(HZ_rect,xhss,1000,500);
% toc
% norm(Z*xmin2-b)
% norm(xmin2)
% norm(xtrue)


%% 

% Ms = [1000,1000,2000,2000,2317,3000,1500,3500];
% Ns = [1500,4000,3000,3500,3141,3600,5000,5000];
Ms = [1000,1500,2000,2500,3000,3500,4000,4500,5000];
Ns = [5001,5001,5001,5001,5001,5001,5001,5001,5001];
k = 50;
blocksize = 200;
samples = 5;
constimings = zeros(size(Ms,2),1);
mltimes = zeros(size(Ms,2),1);
%initsolve_timesn = zeros(size(Ms,2),1);
initsolve_times = zeros(size(Ms,2),1);
savedsolve_times = zeros(size(Ms,2),1);
%mlsolve_errors = zeros(size(Ms,2),1);
hsssolve_errors = zeros(size(Ms,2),1);

for i=1:size(Ms,2)
    M = Ms(i);
    N = Ns(i);
    [Z,HZ,constimings(i)] = rectsampler(M, N, k, blocksize);
    mltime = 0;
    savedsolve_time = 0;
    mlsolve_error = 0;
    hsssolve_error = 0;
    x = rand(N,1);
    b = Z*x;
    % global use_range
    % use_range = 1;
    tic
    x = HZ\b; %#ok<RHSFN>
    initsolve_time = toc;
    % use_range = 0;
    % tic
    % [x,HZ] = HZ\b; %#ok<RHSFN>
    % initsolve_timen = toc;
    for j=1:samples
        x = rand(N,1);
        b = Z*x;
        
        tic
        xtrue = lsqminnorm(Z,b);
        mltime = toc+mltime;
        %mlsolve_error = ml_solve;

        tic
        [x,HZ] = HZ\b; %#ok<RHSFN>
        savedsolve_time = toc + savedsolve_time;
        hsssolve_error = hsssolve_error + norm(x-xtrue)/norm(xtrue);
    end
    constimings
    mltimes(i) = mltime/samples
    % initsolve_timesn(i) = initsolve_timen
    initsolve_times(i) = initsolve_time
    savedsolve_times(i) = savedsolve_time/samples
    hsssolve_errors(i) = hsssolve_error/samples
end

%% Table it
MatrixSize = strcat(string(Ms(:)), " x ", string(Ns(:)));

TotalElements = Ms(:) .* Ns(:);

T = table( ...
    MatrixSize, ...
    TotalElements, ...
    constimings(1:numel(Ms)), ...
    mltimes(1:numel(Ms)), ...
    initsolve_times(1:numel(Ms)), ...
    savedsolve_times(1:numel(Ms)), ...
    hsssolve_errors(1:numel(Ms)), ...
    'VariableNames', ...
    {'Matrix Size','NumElements','HSS Construction Time','Matlab Solve Time','Initial Solve Time (Range)','Stored Solve Time','Relative Solve Error'} ...
);

T.NumElements = T.NumElements / 1e6;
% T.Properties.VariableUnits{2} = 'million';

% T.InitSolveTime = round(1e3*T.InitSolveTime,2);
% T.MLTime        = round(1e3*T.MLTime,2);
% T.SavedSolveTime= round(1e3*T.SavedSolveTime,2);

% T.Properties.VariableUnits = {'','million','ms','ms','ms'};

T = sortrows(T,'NumElements');


disp(T)




%% helper functions
function Z = squaresampler(N)
x_all = linspace(-10,10,2*N);
x = x_all(1:2:2*N);
y = x_all(2:2:2*N);
[X,Y] = meshgrid(x,y);
Z = 1 ./ (X-Y) + eye(N);
end

function [Z,HZ,varargout] = rectsampler(M,N,k,blocksize)
Z = rand(M,k)*rand(k,N);
tic
HZ = hss(Z,blocksize = blocksize);
timing = toc;
if nargout == 3
    varargout{1} = timing;
end
HZ = diagmodify(HZ);
Z = full(HZ);
end

function H = diagmodify(H)
if H.isleaf
    H.D = H.D+rand(size(H.D));
else
    H.A11 = diagmodify(H.A11);
    H.A22 = diagmodify(H.A22);
end
end

function xmin = minnormv1(H,x0,b,iters)
N = size(x0,1);
Null = zeros(N,iters);
for i=1:iters
    v = randn(N,1);
    bk = b+H*v; %nlogn
    xk = H\bk; %
    z = xk-x0-v;
    Null(:,i) = z;
end
% 
% tic
% xmin1 = x0-Null*((Null'*Null)\(Null' * x0));
% toc

[Q,~] = qr(Null,'econ');
xmin = x0-Q*Q'*x0;

end

function xmin = minnormv2(H,x0,size,iters)
M = H.size(1);

Ht = H.';

V = randn(M,size);
G = Ht*V;
[Q,~] = qr(G);
xmin = Q(:,end-iters:end)*(Q(:,end-iters:end))'*x0;
end

function A = Hankel_mat(N)
    t = linspace(0,2*pi,N+1); t(end) = []; z = exp(1i*t(:)); %circle points
    zp = 1i*z;
    w = (2*pi/N)*abs(zp); %weights for nystrom
    d = bsxfun(@minus,z,z.'); % displacements matrix
    r = abs(d);
    ka = 10;
    A = (.25i * besselh(0,ka*r)) .* w;
    diag_ind = sub2ind([N,N],1:N,1:N);
    A(diag_ind) = 1;
end

%% helpers for hankel matrix kernel representation
function A = Hankel_kernel(N)
    % Geometry (stored, but O(N), not O(N^2))
    t = linspace(0, 2*pi, N+1); 
    t(end) = [];
    
    z  = exp(1i*t(:));      % circle points
    zp = 1i*z;
    w  = (2*pi/N)*abs(zp);  % Nystrom weights
    
    ka = 10;
    
    % Kernel handle: returns |I| x |J| matrix
    A = @(I,J) hankel_block(I,J,z,w,ka);
end

function B = hankel_block(I,J,z,w,ka)
    zi = z(I);
    zj = z(J);
    
    % pairwise distances
    d = zi - zj.';      % |I| x |J|
    r = abs(d);
    
    % Helmholtz Hankel kernel
    B = 0.25i * besselh(0, ka*r) .* w(J).';
    
    % Diagonal correction
    [ii,jj] = find(bsxfun(@eq, I(:), J(:).'));
    if ~isempty(ii)
        B(sub2ind(size(B),ii,jj)) = 1;
    end
end

function y = lsqr_wrapper(x, mode, H)
    if strcmp(mode, 'notransp')
        y = H * x; % Ax
    else
        y = H.' * x; % A'x
    end
end

hker = Hankel_kernel(3000);
H = hss(hker);
A = full(H);
cond(A)
