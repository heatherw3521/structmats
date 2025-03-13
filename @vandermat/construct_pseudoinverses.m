function [leftpseudo, rightpseudo, V] = construct_pseudoinverses(V, tol)
%CONSTRUCT_PSEUDOINVERSES Constructs the left and right pseudoinverses of a
% Vandermonde matrix V. The pseudoinverses are computed from a low-rank
% approximation of V within tolerance eps. Note that the ill-conditioning
% of Vandermonde matrices makes inversion a fool's errand for many sets of
% nodes.
arguments
    V;
    tol = 1e-12;
end

% Get basic information
x =  V.x;
n = V.n;

%%% COMPUTE A LOW-RANK APPROXIMATION TO V
if(any(abs(imag(x))>1e-10))
    error("Can't currently handle complex-valued nodes! This feature is coming soon...");
else
    % Get a bound on the numerical rank.
    r = numericalrank_bound(V,tol);
    if(r > max(size(V))) % Nothing to be gained from doing fADI
        [U,s,V] = svd(full(V),'econ','vector');
        leftpseudo = @(B) V * (s .\ (U'*B));
        rightpseudo = @(B) U * (s .\ (V'*B));
    else % Compression is possible!
        % Form the spectral sets to solve Zolotarev 3 on
        E = x; 
        F = exp((1:n)*2i*pi/n).';

        % Solve Zolotarev 3 to get the shift parameters
        [zer, pol] = Z3(E,F,r,r,200,0.95,1e-14);

        % Use fADI to obtain a low-rank approximation
        [Z,d,Y] = fadi(V.A,V.B,V.u,V.v,pol,zer); 
        Z = Z .* d.'; % Absord d into Z so that V ~ ZY'

        % Iterated SVDs to cheaply construct a cheap pseudoinverse
        % Compute SVDs of Z and Y
        [UZ,sZ,VZ] = svd(Z,'econ','vector');
        [UY,sY,VY] = svd(Y,'econ','vector');
        % Compute SVD of the product ZY*
        [UU,s,VV] = svd(sZ .* (VZ'*VY) .* sY','econ','vector');
        U = UZ*UU;
        V = UY*VV;

        % Adhoc stabilization attempt ..?
        % I = s>s(1)*1e-12; 
        % s(I) = 1./s(I);
        % s(~I) = 0;
        % sinv = s;
        
        % ASSEMBLE THE PSEUDO-INVERSES
        leftpseudo = @(B) V*(s .\ (U'*B));
        rightpseudo = @(B) U*(s .\ (V'*B));
    end
    
end

