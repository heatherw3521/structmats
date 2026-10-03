%%
% M = 100;
% N = 200;
% tol = 1e-8;
% x_all = linspace(-10,10,2*N);
% x = x_all(1:2:2*N);
% y = x_all(2:2:2*N);
% [X,Y] = meshgrid(x,y);
% Z = 1 ./ (X-Y) + eye(N);
% 
% 
% rectZ = Z(1:M,1:N);
% 
% rectH = hss(rectZ,blocksize = 50,k=48);
% 
% fullrectH = full(rectH);
% 
% norm(rectZ-fullrectH)/norm(rectZ)

%% 
M = randi(100)+1000;
N = randi(100)+M+1000;
k = randi(50)+50;
blocksize = 200;
X = rand(M,k);
Y = rand(k,N);
rectZ = X*Y;
levelcount = floor(log(max(size(rectZ))/blocksize)/log(2));


A = [eye(floor(M/(2^levelcount))) flip(eye(floor(M/(2^levelcount))))];
%A = padarray(A,[floor(M/(2^levelcount))-size(A,1),floor(N/(2^levelcount))-size(A,2)],0);
Ar = repmat(A, 1, 2^levelcount);                                  
Ac = mat2cell(Ar, size(A,1), repmat(size(A,2),1,2^levelcount));    
Ab = blkdiag(Ac{:});                              
rectZ = rectZ + resize(Ab,size(rectZ));

rectH = hss(rectZ,blocksize = blocksize);


norm(rectZ-full(rectH))/norm(rectZ)

%% 
b = rand(M,1);
tic
xtrue = lsqminnorm(rectZ,b);
toc
%norm(xorig);
%norm(xtrue);
%norm(xtrue-xorig)/norm(xorig);

tic
hssx = rectH\b;
toc

norm(hssx-xtrue)/norm(xtrue)

rectZ*hssx-b;
norm(rectZ*hssx-b)


norm(xtrue)
norm(hssx)
x0 = hssx;

%% lets try to make it min norm solution


%generate null space

Null = [];
count = 15;
for i=1:count
    v = randn(N,1);
    bk = b+rectH*v; %nlogn
    xk = rectH\bk; %
    z = xk-x0-v;
    Null = [Null z];
end

tic
xmin1 = x0-Null*((Null'*Null)\(Null' * x0));
toc

tic
[Q,R] = qr(Null,'econ');
xmin2 = x0-Q*Q'*x0;
toc

norm(xtrue)
norm(hssx)
norm(rectZ*xmin1-b)
norm(xmin1)
norm(rectZ*xmin2-b)
norm(xmin2)


%% round 2 on minnorm solution

% lets make a null space!!
count = max(M,N)-1;

rectHt = rectH.';
% rectZt = rectZ';
% rectHt = hss(rectZt,blocksize = blocksize);

V = randn(M,count);
G = rectHt*V;

% [Q,R] = qr(G,'econ');
[Q,R] = qr(G);
% xmin3 = Q*Q'*x0;
xmin3 = Q(:,end-count:end)*(Q(:,end-count:end))'*x0;

norm(rectZ*xmin3-b)
norm(xtrue)
norm(xmin3)

%%
M = 30;
N = 40;
k = 2;
blocksize = 20;
[Z,HZ] = rectsampler(M,N,k,blocksize);
HZ_orig = HZ;
% one level first

LB = [HZ.A21.Y; HZ.A11.D];
RB = [HZ.A12.Y; HZ.A22.D];

[QL,RL] = qr(LB');
[QR,RR] = qr(RB');

Qs = blkdiag(QL,QR);

HZ.A11.D = RL(1:17,3:end)';
HZ.A21.Y = RL(1:17,1:2)';
HZ.A22.D = RR(1:17,3:end)';
HZ.A12.Y = RR(1:17,1:2)';
fZ = full(HZ);

% Generate a new random vector for the next iteration
b = rand(M, 1);
% Compute the least squares solution with the updated vector
xtrue = lsqminnorm(Z, b);
% Calculate the norm of the difference between the new solution and the previous one
x = lsqminnorm(fZ, b);
x = HZ\b;
x = [x(1:17);zeros(3,1);x(18:end);zeros(3,1)];
x = Qs*x;


%% Helper functions

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