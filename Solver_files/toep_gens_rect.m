function [G, H] = toep_gens_rect(tr, tc, q) 
%given a toeplitz matrix with 
% rowvec t0, colvec t1, 
% this constructs generators GH^*
% so that Z_{q}T - TZ_{1} = GH^*
%%

tr = tr(:); 
tc = tc(:); 
n = length(tr); 
m = length(tc);

if m == n
    error('structmats:toep_gens_rect:use special code for square matrices instead.')
end

tcs = flip(tc(2:end));
r1 = q*tcs(1:n)-[tr(2:end);tr(1)];
H = [r1.';zeros(1,n-1), 1]; 
r2 = [tc(m-n+1); flip(tr);tc(2:m-n)]-tc; 
r2(1)=0;
G = [ [1; zeros(m-1,1)],r2]; 
H = H';
end