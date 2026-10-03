A = randn(40,40);
m = 40;
n = 40;
X = linspace(0,1,n);
Y = linspace(3,4,m).';
C = 1 ./ (X-Y);
A = C(1:m,1:n);
[Z, cols] = inter_decompv2(A, ctype = 'k', cval = 10, orientation= 'columns');
err = norm(A - A(:,cols)*Z, 'fro') / norm(A, 'fro')