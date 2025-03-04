function x = mldivide(T, b)
%MLDIVIDE Left matrix divide for TOEPLITZMAT.
%   Solves the equations Tx = b
%
% functionality is currently limited to square systems with 
% dimensions of powers of 2.
% See Solver_files/structsolv_toeplitz.m for more details
tc = T.tc; 
tr = T.tr;
x = structsolv_toeplitz(tc,tr, b);
%error('Toeplitz left-divide implementation in progress!')
%x = Toeps(T.tc, T.tr, b);
end

