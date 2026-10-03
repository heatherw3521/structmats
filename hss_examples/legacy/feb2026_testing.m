%% Rectangular HSS underdetermined solve scaling benchmark with accuracy
clear; clc;

% Parameters
Ms = round(logspace(3, 4.5, 10));   % Row dimensions
aspect_ratios = [1.5, 2, 2.5, 3.0]; % Different N/M ratios
k = 45;                              % HSS rank parameter
nreps = 3;                          % Number of solve repetitions
save_hss_structs = true;            % Set true to save HSS objects
alpha = 1e-5; % 1e-5 gives a condition numer ~1e10 on cauchy rectsamplerillcond function

% Estimate max size for full matrices (conservative, assuming 8 bytes per double)
max_bytes = 8e9; % 8 GB limit for safety
% max_full_size = floor(sqrt(max_bytes / 8)); % For M*N elements
max_full_size = 15000;

% Storage structure
results = struct();

% Main loop over aspect ratios
for ar_idx = 1:length(aspect_ratios)
    aspect_ratio = aspect_ratios(ar_idx);
    
    fprintf('\n########################################\n');
    fprintf('ASPECT RATIO: %.1f:1\n', aspect_ratio);
    fprintf('########################################\n');
    
    N_sizes = round(Ms * aspect_ratio);
    
    % Initialize storage for this aspect ratio
    hss_build_time = zeros(size(Ms));
    matlab_lsqminnorm_time = zeros(size(Ms));
    hss_backslash_time = zeros(size(Ms));
    hss_factored_time = zeros(size(Ms));
    hss_lsqr_time = zeros(size(Ms));
    
    % Accuracy metrics (2-norm of residual)
    matlab_residual_norm = zeros(size(Ms));
    hss_backslash_residual_norm = zeros(size(Ms));
    hss_factored_residual_norm = zeros(size(Ms));
    hss_lsqr_residual_norm = zeros(size(Ms));
    
    % Storage metrics (in MB)
    full_storage = zeros(size(Ms));
    hss_storage = zeros(size(Ms));
    hss_factored_storage = zeros(size(Ms));
    
    if save_hss_structs
        HSS_cells = cell(size(Ms));
    end
    
    % Loop over sizes
    for idx = 1:length(Ms)
        M = Ms(idx);
        N = N_sizes(idx);
        
        fprintf('\n====================================\n');
        fprintf('Size: M = %d, N = %d (%.1f:1)\n', M, N, N/M);
        fprintf('====================================\n');
        
        % Check if we can form full matrix
        can_form_full = (M * N < max_full_size^2);
        % --- HSS Construction ---
        blocksize = min(max(200, round(0.5*sqrt(M))), 2000);
        fprintf('Building HSS with blocksize %d...\n', blocksize);
        %[Z_handle, HZ_rect, build_time] = rectsamplerillcond(M, N, k, blocksize,'condition_number', condnum);
        %[Z_handle, HZ_rect, build_time] = rectsamplerillcond(M, N, k, blocksize, 'alpha', alpha, 'form_full', can_form_full);
        [Z_handle, HZ_rect, build_time] = rectsampler(M, N, k, blocksize);
        if can_form_full
            Z = full(HZ_rect);
            %condnum = cond(Z);
        end

        hss_build_time(idx) = build_time;
        fprintf('HSS build time: %.3f s\n', hss_build_time(idx));
        
        % Measure HSS storage
        hss_info = whos('HZ_rect');
        hss_storage(idx) = hss_info.bytes / 1e6; % Convert to MB
        
        if save_hss_structs
            HSS_cells{idx} = HZ_rect;
        end
        
        % Create test problem with known solution
        x_true = randn(N, 1);
        x_true = x_true / norm(x_true); % Normalize
        if can_form_full
            b = Z * x_true;
        else
            b = HZ_rect *x_true;
        end
        
        % --- Method 1: MATLAB lsqminnorm (only if full matrix fits) ---
        if can_form_full
            fprintf('Testing MATLAB lsqminnorm...\n');
            
            % Measure full matrix storage
            full_info = whos('Z');
            full_storage(idx) = full_info.bytes / 1e6; % MB
            
            tic;
            for r = 1:nreps
                %x1 = lsqminnorm(Z, b);
                x1 = lsqr(Z,b,0,1000);
            end
            matlab_lsqminnorm_time(idx) = toc / nreps;
            
            % Compute accuracy (2-norm of residual)
            matlab_residual_norm(idx) = norm(Z*x1 - b);
            
            fprintf('Avg MATLAB lsqminnorm: %.6f s | ||Ax-b||: %.2e\n', ...
                    matlab_lsqminnorm_time(idx), matlab_residual_norm(idx));
            
            % Clear full matrix after MATLAB test
            clear Z
        else
            fprintf('Skipping MATLAB lsqminnorm (matrix too large: %d x %d)\n', M, N);
            matlab_lsqminnorm_time(idx) = NaN;
            matlab_residual_norm(idx) = NaN;
            
            % Estimate full storage
            full_storage(idx) = (M * N * 8) / 1e6; % MB
        end
        
        % % --- Method 2: HSS backslash (unfactored) ---
        % %can_form_full = false;
        % if can_form_full
        %     fprintf('Testing HSS backslash...\n');
        %     tic;
        %     for r = 1:nreps
        %         [x2, HZ_factored] = HZ_rect \ b;
        %     end
        %     hss_backslash_time(idx) = toc / nreps;
        % 
        %     % Compute accuracy using HSS matvec
        %     hss_backslash_residual_norm(idx) = norm(HZ_rect*x2 - b);
        % 
        %     fprintf('Avg HSS backslash: %.6f s | ||Ax-b||: %.2e\n', ...
        %         hss_backslash_time(idx), hss_backslash_residual_norm(idx));
        % 
        %     % Measure factored HSS storage
        %     factored_info = whos('HZ_factored');
        %     hss_factored_storage(idx) = factored_info.bytes / 1e6; % MB
        % 
        %     % Check if factored HSS fits in memory
        %     can_use_factored = (hss_factored_storage(idx) * 1e6 < max_bytes);
        % else
        %     fprintf('Skipping HSS backslash (matrix too large: %d x %d)\n', M, N);
        %     hss_backslash_time(idx) = NaN;
        %     hss_backslash_residual_norm(idx) = NaN;
        %     can_use_factored = 0;
        % end
        % 
        % if can_use_factored
        %     % --- Method 3: HSS with stored factors ---
        %     fprintf('Testing HSS with stored factors...\n');
        %     tic;
        %     for r = 1:nreps
        %         x3 = HZ_factored \  b;
        %     end
        %     hss_factored_time(idx) = toc / nreps;
        % 
        %     % Compute accuracy using HSS matvec
        %     hss_factored_residual_norm(idx) = norm(HZ_rect*x3 - b);
        % 
        %     fprintf('Avg HSS factored: %.6f s | ||Ax-b||: %.2e\n', ...
        %             hss_factored_time(idx), hss_factored_residual_norm(idx));
        % else
        %     fprintf('Skipping HSS factored solve (factored storage %.1f MB exceeds limit)\n', ...
        %             hss_factored_storage(idx));
        %     hss_factored_time(idx) = NaN;
        %     hss_factored_residual_norm(idx) = NaN;
        %     hss_factored_storage(idx) = full_storage(idx);
        %     % Clear factored HSS to free memory
        %     clear HZ_factored
        % end
        
        % --- Method 4: HSS LSQR ---
        fprintf('Testing HSS LSQR...\n');
        fun = @(x, mode) lsqr_wrapper(x, mode, HZ_rect);
        tic;
        for r = 1:nreps
            [x4, flag] = lsqr(fun, b, 0, 1000);
        end
        hss_lsqr_time(idx) = toc / nreps;
        
        % Compute accuracy using HSS matvec
        hss_lsqr_residual_norm(idx) = norm(HZ_rect*x4 - b);
        
        fprintf('Avg HSS LSQR: %.6f s | ||Ax-b||: %.2e\n', ...
                hss_lsqr_time(idx), hss_lsqr_residual_norm(idx));
        
        % Print storage summary
        fprintf('\nStorage: Full=%.1f MB | HSS=%.1f MB | HSS+factors=%.1f MB\n', ...
                full_storage(idx), hss_storage(idx), hss_factored_storage(idx));
        if can_form_full && can_use_factored
            fprintf('Compression ratio: Full/HSS = %.1fx | Full/HSS+factors = %.1fx\n', ...
                    full_storage(idx)/hss_storage(idx), ...
                    full_storage(idx)/hss_factored_storage(idx));
        elseif can_use_factored
            fprintf('Compression ratio: Full/HSS = %.1fx (estimated) | Full/HSS+factors = %.1fx\n', ...
                    full_storage(idx)/hss_storage(idx), ...
                    full_storage(idx)/hss_factored_storage(idx));
        end
        
        % Clean up large objects
        clear HZ_rect HZ_factored x2 x3 x4
    end
    
    % Store results for this aspect ratio
    results(ar_idx).aspect_ratio = aspect_ratio;
    results(ar_idx).M = Ms(:);
    results(ar_idx).N = N_sizes(:);
    results(ar_idx).hss_build_time = hss_build_time(:);
    results(ar_idx).matlab_time = matlab_lsqminnorm_time(:);
    results(ar_idx).hss_backslash_time = hss_backslash_time(:);
    results(ar_idx).hss_factored_time = hss_factored_time(:);
    results(ar_idx).hss_lsqr_time = hss_lsqr_time(:);
    results(ar_idx).matlab_residual = matlab_residual_norm(:);
    results(ar_idx).hss_backslash_residual = hss_backslash_residual_norm(:);
    results(ar_idx).hss_factored_residual = hss_factored_residual_norm(:);
    results(ar_idx).hss_lsqr_residual = hss_lsqr_residual_norm(:);
    results(ar_idx).full_storage = full_storage(:);
    results(ar_idx).hss_storage = hss_storage(:);
    results(ar_idx).hss_factored_storage = hss_factored_storage(:);
    
    if save_hss_structs
        results(ar_idx).HSS_cells = HSS_cells;
    end
    
    % Create table for this aspect ratio
    results(ar_idx).table = table(Ms(:), N_sizes(:), ...
                                  hss_build_time(:), ...
                                  matlab_lsqminnorm_time(:), ...
                                  hss_backslash_time(:), ...
                                  hss_factored_time(:), ...
                                  hss_lsqr_time(:), ...
                                  matlab_residual_norm(:), ...
                                  hss_backslash_residual_norm(:), ...
                                  hss_factored_residual_norm(:), ...
                                  hss_lsqr_residual_norm(:), ...
                                  full_storage(:), ...
                                  hss_storage(:), ...
                                  hss_factored_storage(:), ...
        'VariableNames', {'M', 'N', 'HSS_build', ...
                          'MATLAB_time', 'HSS_backslash_time', 'HSS_factored_time', 'HSS_lsqr_time', ...
                          'MATLAB_residual', 'HSS_backslash_residual', 'HSS_factored_residual', 'HSS_lsqr_residual', ...
                          'Full_storage_MB', 'HSS_storage_MB', 'HSS_factored_storage_MB'});
    
    % Save individual table
    filename = sprintf('hss_rect_solve_AR%.1f.txt', aspect_ratio);
    writetable(results(ar_idx).table, filename, 'Delimiter', '\t');
end

% Save all results
save('hss_rect_solve_results_illcondv3.mat', 'results', 'aspect_ratios');
fprintf('\nBenchmark complete. Results saved.\n');

%% Plot timing results - separate plots for each aspect ratio
load hss_rect_solve_results_illcond.mat

num_ratios = length(aspect_ratios);

for ar_idx = 1:num_ratios
    figure('Position', [100 + 50*ar_idx, 100 + 50*ar_idx, 1000, 700]);

    
    M = results(ar_idx).M;
    
    % MATLAB lsqminnorm
    valid_matlab = ~isnan(results(ar_idx).matlab_time);
    if any(valid_matlab)
        loglog(M(valid_matlab), results(ar_idx).matlab_time(valid_matlab), ...
               'bo-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
               'DisplayName', 'MATLAB lsqminnorm');
    end
    hold on;
    % HSS Backslash
    % loglog(M, results(ar_idx).hss_backslash_time, ...
    %        's-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
    %        'DisplayName', 'HSS Backslash');

    % HSS Factored
    % valid_factored = ~isnan(results(ar_idx).hss_factored_time);
    % if any(valid_factored)
    %     loglog(M(valid_factored), results(ar_idx).hss_factored_time(valid_factored), ...
    %            'd-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
    %            'DisplayName', 'HSS Factored');
    % end
    
    colors = lines(4);
    purple = colors(4,:);
    % HSS LSQR
    loglog(M, results(ar_idx).hss_lsqr_time, ...
           '^-', 'Color', purple, 'LineWidth', 2.5, 'MarkerSize', 8, ...
           'DisplayName', 'HSS LSQR');
    % Add reference curve
    mid_idx = ceil(length(M)/2);
    ref_nlogn = M .* log(M);
    c = results(ar_idx).hss_lsqr_time(mid_idx) / ref_nlogn(mid_idx);
    loglog(M, c * ref_nlogn, 'k--', 'LineWidth', 2, 'DisplayName', 'O(M log M)');
    
    grid on;
    xlabel('M (number of rows)', 'FontSize', 14);
    ylabel('Solve Time (seconds)', 'FontSize', 14);
    title(sprintf('Underdetermined Solve: Aspect Ratio %.1f:1 (N = %.1f×M)', ...
                  aspect_ratios(ar_idx), aspect_ratios(ar_idx)), 'FontSize', 16);
    legend('Location', 'NorthWest', 'FontSize', 12);
    set(gca, 'FontSize', 12);
    
    hold off;
end

%% Plot residual error results - separate plots for each aspect ratio
load hss_rect_solve_resultsv2.mat

% Only use first 2 aspect ratios
num_ratios = min(2, length(aspect_ratios));

for ar_idx = 1:num_ratios
    figure('Position', [100 + 50*ar_idx, 100 + 50*ar_idx, 1000, 700]);
    hold on;
    
    M = results(ar_idx).M;
    
    % MATLAB lsqminnorm
    valid_matlab = ~isnan(results(ar_idx).matlab_residual);
    if any(valid_matlab)
        loglog(M(valid_matlab), results(ar_idx).matlab_residual(valid_matlab), ...
               'o-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
               'DisplayName', 'MATLAB lsqminnorm');
    end
    
    % HSS Backslash
    loglog(M, results(ar_idx).hss_backslash_residual, ...
           's-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
           'DisplayName', 'HSS Backslash');
    
    % HSS Factored
    valid_factored = ~isnan(results(ar_idx).hss_factored_residual);
    if any(valid_factored)
        loglog(M(valid_factored), results(ar_idx).hss_factored_residual(valid_factored), ...
               'd-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
               'DisplayName', 'HSS Factored');
    end
    
    % HSS LSQR
    loglog(M, results(ar_idx).hss_lsqr_residual, ...
           '^-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
           'DisplayName', 'HSS LSQR');
    
    grid on;
    xlabel('M (number of rows)', 'FontSize', 14);
    ylabel('Residual Norm ||Ax-b||_2', 'FontSize', 14);
    title(sprintf('Residual Error: Aspect Ratio %.1f:1 (N = %.1f×M)', ...
                  aspect_ratios(ar_idx), aspect_ratios(ar_idx)), 'FontSize', 16);
    legend('Location', 'Best', 'FontSize', 12);
    set(gca, 'FontSize', 12);
    
    hold off;
end


%% Plot storage results - separate plots for each aspect ratio
load hss_rect_solve_resultsv2.mat

% Only use first 2 aspect ratios
num_ratios = min(2, length(aspect_ratios));

for ar_idx = 1:num_ratios
    figure('Position', [100 + 50*ar_idx, 100 + 50*ar_idx, 1000, 700]);
    hold on;
    
    M = results(ar_idx).M;
    
    % Full matrix storage
    loglog(M, results(ar_idx).full_storage, ...
           'o-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
           'DisplayName', 'Full Matrix');
    
    % HSS storage
    loglog(M, results(ar_idx).hss_storage, ...
           's-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
           'DisplayName', 'HSS');
    
    % HSS + Factors storage
    loglog(M, results(ar_idx).hss_factored_storage, ...
           'd-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
           'DisplayName', 'HSS + Factors');
    
    grid on;
    xlabel('M (number of rows)', 'FontSize', 14);
    ylabel('Storage (MB)', 'FontSize', 14);
    title(sprintf('Storage Requirements: Aspect Ratio %.1f:1 (N = %.1f×M)', ...
                  aspect_ratios(ar_idx), aspect_ratios(ar_idx)), 'FontSize', 16);
    legend('Location', 'NorthWest', 'FontSize', 12);
    set(gca, 'FontSize', 12);
    
    hold off;
end


%% spy stuff

M = 200;
N = M*2;
k = 50;
blocksize = 200;
[Z,HZ_rect,time] = rectsampler(M, N, k, blocksize);
disp("----------")
disp("Size = " + M + " x " + N)
disp("HSS construction: " + time + " seconds")
figure
spy(HZ_rect)

%% MAT MAT TESTING (small first)
m1 = 20;
n = 20;
m2 = 20;

A1 = randn(m1,n);
A2 = randn(n,m2);

H1 = hss(A1,blocksize = 10);
H2 = hss(A2,blocksize = 10);

A = A1*A2;
% H = H1*H2;

H = hss_matmatv2(H1,H2);



%% Helper functions
function [Z_handle, HZ, varargout] = rectsampler(M, N, k, blocksize, varargin)
    % Create structured low-rank rectangular matrix
    % Optional: kernel_type = 'random' (default), 'oscillatory', 'decay'
    
    p = inputParser;
    addParameter(p, 'kernel_type', 'random', @ischar);
    addParameter(p, 'freq', 10, @isnumeric);  % for oscillatory
    addParameter(p, 'decay_rate', 0.1, @isnumeric);  % for decay
    parse(p, varargin{:});
    
    kernel_type = p.Results.kernel_type;
    
    % Create geometry (like your Hankel example)
    t_row = linspace(0, 1, M+1); t_row(end) = [];
    t_col = linspace(0, 1, N+1); t_col(end) = [];
    
    % Store geometry for block function
    geom.t_row = t_row(:);
    geom.t_col = t_col(:);
    geom.k = k;
    geom.kernel_type = kernel_type;
    geom.freq = p.Results.freq;
    geom.decay_rate = p.Results.decay_rate;
    
    % Create function handle
    Z_handle = @(I, J) rect_kernel_block(I, J, geom);
    
    % Build HSS from function handle
    tic;
    HZ = hss(Z_handle, blocksize=blocksize, sizeA=[M, N]);
    timing = toc;
    
    % Apply diagonal modification
    HZ = diagmodify(HZ);
    
    % Return timing if requested
    if nargout == 3
        varargout{1} = timing;
    end
end

function [Z, HZ, varargout] = rectsamplerillcond(M, N, k, blocksize, varargin)
    p = inputParser;
    addParameter(p, 'alpha', 1e-5, @isnumeric);
    addParameter(p, 'form_full', 0);
    parse(p, varargin{:});
    
    alpha = p.Results.alpha;
    form_full = p.Results.form_full;
    
    s = linspace(0, 1, M)';
    t = linspace(0, 1, N)';
    
    Z_handle = @(I,J) exp(-alpha * abs(s(I) - t(J).'));
    
    tic;
    HZ = hss(Z_handle, blocksize=blocksize, sizeA=[M, N]);
    timing = toc;
    
    if form_full
        Z = full(HZ);
    else
        Z = [];
    end
    
    if nargout == 3
        varargout{1} = timing;
    end
end

function B = rect_kernel_block(I, J, geom)
    % Compute block of low-rank structured matrix
    
    ti = geom.t_row(I);
    tj = geom.t_col(J);
    k = geom.k;
    
    switch geom.kernel_type
        case 'random'
            % Random low-rank (fastest, simplest)
            % Generate deterministically based on indices for consistency
            rng(12345);  % Fixed seed for reproducibility
            U_full = randn(length(geom.t_row), k) / sqrt(k);
            V_full = randn(length(geom.t_col), k) / sqrt(k);
            B = U_full(I, :) * V_full(J, :)';
            
        case 'oscillatory'
            % Oscillatory kernel (like discretized integral operator)
            % B(i,j) = sum_{l=1}^k sin(freq*l*ti) * cos(freq*l*tj)
            freq = geom.freq;
            U = zeros(length(I), k);
            V = zeros(length(J), k);
            for l = 1:k
                U(:, l) = sin(freq * l * ti) / sqrt(k);
                V(:, l) = cos(freq * l * tj) / sqrt(k);
            end
            B = U * V';
            
        case 'decay'
            % Exponentially decaying kernel
            % Mimics distance-based decay
            [Ti, Tj] = meshgrid(tj, ti);
            dist = abs(Ti - Tj);
            
            % Low-rank approximation via truncated SVD of decay kernel
            % For small blocks, just compute directly
            if length(I)*length(J) < 1000
                B = exp(-geom.decay_rate * dist);
            else
                % For larger blocks, use low-rank structure
                U = zeros(length(I), k);
                V = zeros(length(J), k);
                for l = 1:k
                    lambda = geom.decay_rate * l;
                    U(:, l) = exp(-lambda * ti) / sqrt(k);
                    V(:, l) = exp(-lambda * tj) / sqrt(k);
                end
                B = U * V';
            end
            
        case 'polynomial'
            % Polynomial kernel: (1 + ti*tj')^k approximation
            U = zeros(length(I), k);
            V = zeros(length(J), k);
            for l = 1:k
                U(:, l) = ti.^(l-1) / sqrt(factorial(l-1));
                V(:, l) = tj.^(l-1) / sqrt(factorial(l-1));
            end
            B = U * V';
            
        otherwise
            error('Unknown kernel type: %s', geom.kernel_type);
    end
end

function H = diagmodify(H)
if H.isleaf
    %H.D = H.D+rand(size(H.D));
    [Q,~] = qr(randn(size(H.D)));
    dg = 20.^(-0:size(H.D,1)-1);
    blk = Q*diag(dg)*Q';
    newblk = zeros(size(H.D));
    newblk(1:size(blk,1),1:size(blk,2)) = blk;
    H.D = H.D + newblk;
else
    H.A11 = diagmodify(H.A11);
    H.A22 = diagmodify(H.A22);
end
end

function y = lsqr_wrapper(x, mode, H)
    if strcmp(mode, 'notransp')
        y = H * x; % Ax
    else
        y = H.' * x; % A'x
    end
end

function HZ = rectsamplerillcondv2(M,N,k,blocksize)
dg = 20.^(-0:10);
m = 10; 
dg = 20.^(-0:m-1);
D = diag(dg);
diagblks = {}; 
for k = 1:6
    [Q, ~] = qr(randn(m,m));
    diagblks{k} = Q*D*Q'; 
end
MD = blkdiag(diagblks{:});
[a, b] = size(MD);
k = 3;
% 
MB = randn(a,k)*randn(k,a);
M = MD + MB; 
M = M(1:4*m,:);
[a, b] = size(M);
%
s = svd(M);
%
B = randn(a,1);
x = lsqr(M,B, 0, 1000);
end