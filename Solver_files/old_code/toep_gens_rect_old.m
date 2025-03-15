function [G, H] = toep_gens_rect_old(tr, tc) 
%given a toeplitz matrix with 
% rowvec t0, colvec t1, 
% this constructs generators GH^*
% so that Q_{1}T - TQ_{-1} = GH^*
%%

tr = tr(:); 
tc = tc(:); 
n = length(tr); 
m = length(tc);
tcs = flip(tc(2:end));
r1 = [tcs(1:n-1);0]-[tr(2:end);0];
%r1 = [-tcs(1:n-1);0]-[tr(2:end);0];
H = [zeros(1,n-1), 1; r1.'];

%r2 = flip(tr(2:end)) + tc(2:end);  
if m==n
    %r2 = flip(flip(tc)- [tr(2:end); -tr(1)]);
    r2 = [2*tr(1); flip(tr(2:end))] + [0;tc(2:end)];
else
    r2 = [tc(m-n+1); flip(tr);tc(2:m-n)]+tc; 
end
%r2 = [0; r2]; 
%eye = [1; zeros(n-1, 1)];
%G = [ eye r2]; 
G = [r2, [1; zeros(m-1,1)]]; 
%H = conj(H); % so that GH^* gets the right answer
H = H.';
end