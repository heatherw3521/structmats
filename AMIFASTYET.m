% Levent Batacki 
%
% I wonder if any of this is even fast... with all the overhead


%% GENERATE TEST MATRICES
k = 1;
N = 20000:1000:24000;

T = {}; Tm = {}; X = {};
i=1;
for n = N
    T{i} = SampleToeplitz(n,n);
    Tm{i} = full(T{i});
    X{i} = rand(n,k);
    i=i+1;
end

%% RUN TIMING TEST
i=1; T1 = []; T2 = [];
for n = N

    t1 = tic();
    y1 = T{i} * X{i};
    tf1 = toc(t1);    

    t2 = tic();
    y2 = Tm{i} * X{i};
    tf2 = toc(t2);


    T1(i) = tf1;
    T2(i) = tf2;
    disp(n)
   % disp(norm(y1-y2))
    i=i+1;
end

figure(1)
clf
semilogy(N,T1)
hold on;
semilogy(N,T2)
legend("toeplitzmat", "matlab")
disp("MAX BENEFIT: " + T2(end)/T1(end))