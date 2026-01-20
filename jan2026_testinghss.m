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