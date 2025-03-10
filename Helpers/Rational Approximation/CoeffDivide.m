function [C,Q] = CoeffDivide(c,z,k,f)
%COEFFCONDITIONS Gets the top k+1 coefficients of p(x)/(x-z), where p(x) 
%   is given by p(x) = c0 + c1x + .. + cn x^n.
%   
arguments
    c (:,1) {mustBeNumeric}
    z (:,1) {mustBeNumeric}

    % optional
    k = length(z)-1;
    f = ones([length(z) 1])
end

% Make sure we're working with vectors
z = reshape(z,[],1);
c = reshape(c,[],1);
f = reshape(f,[],1);

C = ones(size(z));
n = length(c);
for d =  n-1:-1:n-k
    cc = c(d)-z.*C(:,end);
    C = [C, cc];
    % vv = cc.*f.';
    % for j = 1:size(Q,2)
    %     vv = vv - Q*(Q.'*vv);
    % end
    % Q = [Q, vv/norm(vv)];
end

[Q,~] = qr(C.*f,'econ');