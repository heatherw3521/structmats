%% 1. Lets create a couple matrices that accept HSS structure

%% a) Square 
N = 1000;
Z = squaresampler(N);
x = rand(N,1);

% create the HSS representation
tic
HZ = hss(Z)
toc

% test a matvec
tic
Z*x;
toc

tic
HZ*x;
toc

norm(Z*x-HZ*x)

%% lets up the size of the square
N = 10000;
Z = squaresampler(N);
x = rand(N,1);
tic
HZ = hss(Z);
toc

% test the matvec
tic
Z*x;
toc

tic
HZ*x;
toc

norm(Z*x-HZ*x)


%% well that took a while
% lets be more precise with creating HZ
tic
HZ = hss(Z,blocksize = 500);
toc

% test the matvec
tic
Z*x;
toc

tic
HZ*x;
toc

norm(Z*x-HZ*x)


%% Next, the ULV solver!
N = 5000;
Z5k = squaresampler(N);
b = rand(N,1);
tic
xtrue = Z5k\b;
toc

tic
HZ5k = hss(Z5k,blocksize = 500,tol=1e-10);
toc

tic
xhss = HZ5k\b;
toc

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
tic
xtrue = Z\b;
toc

tic
xhss = HZ\b;
toc

norm(xhss-xtrue)/norm(xtrue)


%% Squares seem to be doing ok. I'm a fan of rectangles though, what gives?

% M = randi(1000)+1000;
% N = randi(1000)+M+3000;
M = 700;
N = 1500;
k = 50;
blocksize = 200;
%Z = rectsampler(M, N, k, blocksize);
Z = rand(M,k)*rand(k,N);
tic
HZ_rect = hss(Z,blocksize = blocksize);
toc
HZ_rect = diagmodify(HZ_rect);
Z = full(HZ_rect);

% tic
% HZ_rect = hss(Z,blocksize=blocksize);
% toc
%%
x = rand(N,1);
b = Z*x;

tic
xtrue = lsqminnorm(Z,b);
toc


tic
[xhss,HZ_rect] = HZ_rect\b; %#ok<RHSFN>
%xhss = minnormv1(HZ_rect,xhss,b,500);
toc

norm(Z*xhss-b)
norm(Z*xtrue-b)
norm(xhss)
norm(xtrue)

%% lets min norm it
x0 = xhss;
H = HZ_rect;
iters = 500;

timing = zeros(iters,1);
error = zeros(iters,1);
xnorm = zeros(iters,1);

% for i=1:iters
%     tic
%     xmin1 = minnormv1(HZ_rect,xhss,b,i);
%     timing(i) = toc;
%     error(i) = norm(Z*xmin1-b);
%     xnorm(i) = norm(xmin1);
% end
N = size(x0,1);
Null = zeros(N,iters);
for i=1:iters
    v = randn(N,1);
    bk = b+H*v; %nlogn
    xk = H\bk; %
    z = xk-x0-v;
    Null(:,i) = z;
end

for i=1:iters
    [Q,~] = qr(Null(:,1:i),'econ');
    xmin = x0-Q*Q'*x0;
    if norm(xmin-xtrue)<1e-4
        closeidx = i;
    end
    xnorm(i) = norm(xmin);
end


%% plot timing, error and xnorm
% plot(1:iters,timing)

plot(1:iters,xnorm,'-',LineWidth=1)
hold on
plot(1:iters,repelem(norm(xtrue),iters),'--',Color='r',LineWidth=0.1)
xline(closeidx)
xlabel('Iterations')
ylabel('norm')
legend('norm(x)','norm(xtrue)')

% tic
% xmin2 = minnormv2(HZ_rect,xhss,1000,500);
% toc
% norm(Z*xmin2-b)
% norm(xmin2)
% norm(xtrue)


%%

% Ms = [1000,1000,2000,2000,2317,3000,1500,3500];
% Ns = [1500,4000,3000,3500,3141,3600,5000,5000];
Ms = [1000,1500,2000,2500,3000];
Ns = [3001,3001,3001,3001,3001];
k = 50;
blocksize = 200;
samples = 5;
constimings = zeros(size(Ms,2),1);
mltimes = zeros(size(Ms,2),1);
initsolve_timesn = zeros(size(Ms,2),1);
initsolve_timesr = zeros(size(Ms,2),1);
savedsolve_times = zeros(size(Ms,2),1);

for i=1:size(Ms,2)
    M = Ms(i);
    N = Ns(i);
    [Z,HZ,constimings(i)] = rectsampler(M, N, k, blocksize);
    mltime = 0;
    savedsolve_time = 0;
    x = rand(N,1);
    b = Z*x;
    global use_range
    use_range = 1;
    tic
    x = HZ\b; %#ok<RHSFN>
    initsolve_timer = toc;
    use_range = 0;
    tic
    [x,HZ] = HZ\b; %#ok<RHSFN>
    initsolve_timen = toc;
    for j=1:samples
        x = rand(N,1);
        b = Z*x;
        
        tic
        xtrue = Z\b;
        mltime = toc+mltime;
    
        tic
        [x,HZ] = HZ\b; %#ok<RHSFN>
        savedsolve_time = toc + savedsolve_time;
    end
    constimings
    mltimes(i) = mltime/samples
    initsolve_timesn(i) = initsolve_timen/samples
    initsolve_timesr(i) = initsolve_timer/samples
    savedsolve_times(i) = savedsolve_time/samples
end

%% Table it
MatrixSize = strcat(string(Ms(:)), " x ", string(Ns(:)));

TotalElements = Ms(:) .* Ns(:);

T = table( ...
    MatrixSize, ...
    TotalElements, ...
    constimings(1:numel(Ms)), ...
    mltimes(1:numel(Ms)), ...
    initsolve_timesn(1:numel(Ms)), ...
    initsolve_timesr(1:numel(Ms)), ...
    savedsolve_times(1:numel(Ms)), ...
    'VariableNames', ...
    {'Matrix Size','NumElements','HSS Construction Time','Matlab Solve Time','Initial Solve Time (Null Space)','Initial Solve Time (Range)','Stored Solve Time'} ...
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
