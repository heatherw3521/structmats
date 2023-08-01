function [THSS, errHSS, cnd] = driver_testinudft(m,n, perturb, plot_pts, get_COND)
% tests NUDFT inversion code on an m by n NUDFT matrix of 
% type-II. For now, choose m > n, m/n = integer. 
%
% V(j,k) = gamma_{j}^{k-1}, gamma = exp(2*pi*1i*p_j/m),
% perturb = max possible perturbation: | p_l -j| < perturb, 
% where l = (j-1)*m/n+1: j*m/n
%plot_pts = 1, plots the nodes. plot_pts = 0, no plot. 
% get_COND = 1, return condition # for V. get_COND = 0 returns []

%%
%ru = linspace(0, 1, n+1); ru = ru(2:end); % start with nth  roots of unity.
%ru = exp(-2*pi*1i*ru);
v = m/n; %for now should be an integer. 

% construct nodes: 
for j = 1:n
    pk = -perturb + 2*perturb*rand(v,1);
    p = j + pk; 
    gamma((j-1)*v+1: j*v) = exp(-2*pi*1i*(p+j-1)/n);
end

%%
if plot_pts
    plot(gamma, '.', 'markersize', 5); 
    hold on
    plot(exp(1i*linspace(0,2*pi, 500)), 'k-')
    hold off
end
%%
% get rhs: 
b = rand(m, 1); 
% get true soln: 
gamma = gamma.'; 
V = gamma.^(0:n-1); 
x = V\b; 
if get_COND
    cnd = cond(V); %get condition number of V. 
else
    cnd = []; 
end
%%
% solve linear system Vx = b using HSS: 
T = tic;
xHSS = INUDFT(gamma,n,b);
THSS = toc(T); 
errHSS = max(abs(x-xHSS));
end