function [r,pol,zer,z,numerator,denominator,f,w,errvec] = RationalApprox(Z,F,m,n,sign_flag,forZ3,nlawson,damping,tol,zolotol)
%RationalApprox Obtains a type (m,n) rational approximation to a given
%   function on a given grid. This is pretty much just the aaa algorithm,
%   written my way! 
% 
%   Some of the code is directly copy-pasted from aaa.
%   My only real contribution is treating the m~=n case.
%   Credit to aaa: https://www.chebfun.org/examples/approx/AAAApprox.html
arguments
    Z,
    F
    
    m = length(Z);
    n = m;
    sign_flag = false;
    forZ3 = false;
    nlawson = 200;
    damping = 0.95;
    tol = 1e-13;
    zolotol = -1;
end
    % Fix the off-by-one error
    mm=m+1;nn=n+1; % pesky thing

    % Ensure we're working with column vectors.
    Z = reshape(Z,[],1);
    F = reshape(F,[],1);

    % ALGORITHM OUTLINE
    % Let's abstract as much as possible. The steps of this algorithm are
    % as follows.
    %
    % We iteratively build a type (m,n) rational function r(Z) in 
    % barycentric form. While type (m,n) is not yet reached,
    %
    % 1. Compute error of approximation over grid
    % 2. Select, greedily, the next 'support' point. We will interpolate at
    %       this point.
    % 3. recalculate the barycentric weights to minimize the error. This
    %     amounts to solving a constrained optimization problem. The
    %     constraints come from the discrepancy between m and n.

    % Ensure the points are unique
    [Z,J] = unique(Z); F = F(J);

    % Initialize weights, supports, function vals, and r(z)
    w = []; z = []; f = [];
    r = @(z) mean(F);

    % Get the discrepancy between m and n
    D = abs(mm-nn);
    c=0;
    if(D > 0) % we will have to constrain the weights
        % c stores the top D coefficients of the polynomial p(z) given by 
        % p(z) = (z-z1)...(z-zk), where z1,...,zk are the support points
        % at iteration k. 
        c = [1 zeros([1 D])]'; 
    end 
    % The weights may be constrained to live in a certain manifold. The
    % matrix C's range is that space's orthogonal complement. The idea is
    % then to use CC to construct an orthogonal projector on the desired
    % space.
    CC = [];

    % Initialize A, which is a big part of the optimization problem. In
    % effect, we solve problems of the form min|Aw|/|w|, with v sometimes
    % being constrained in a certain space
    A = [];

    % Make a copy of the grid and function values. These will be trimmed
    % every time we choose the next support point.
    grid = Z;
    if(zolotol > 0)
        zoloE = Z(Z==-1);
        zoloF = Z(Z==1);
    end
    fvals = F;

    % Check convergence
    abstol = tol * norm(F,inf);
    errvec = [];

    % Build type (k,k), using iterate. Then, constrain to ensure neither of
    % the degrees are too high. 
    % Q: Would it be OK to just constrain at the end?
    for k = 1:max(mm,nn)
        % Choose the next support point
        [zk,fk,j] = select_next_support(grid,fvals,r);
        % Update list of supports and interp. values
        z = [zk;z]; f = [fk;f];
        % Remove the chosen support from future consideration
        grid(j) = []; fvals(j) = [];

        % Update top D polynomial coefficients
        if(D > 0)
           c = update_polynomial_coeffs(c,zk,k,D);
        end

        % Update the weights by solving an optimization problem
        A = (f.'- fvals) ./ (grid - z.');
        [Q,P] = get_constraint(c,k,mm,nn,f,z,forZ3); % Determine constraint
        w = new_weights(A, P, Q, sign_flag); % Solve for new weights
        
        % Compute updated rational approximation
        numerator = @(x) sum((f.*w).' ./ (x-z.'),2);
        denominator = @(x) sum(w.' ./(x-z.'),2);
        r = @(x)  numerator(x) ./ denominator(x);
        
        % Convergence check
        maxerr = norm(fvals - r(grid), inf);
        errvec = [errvec; maxerr];
        if(zolotol<0 && maxerr <= abstol)
            break
        elseif( maxerr<= 1 && (maxerr/(1+sqrt(1-maxerr^2)))^2 < zolotol)
            break
        end
    end

    % Lawson iteration to get new weights
    if(nlawson>0)
        [ww, maxerr_lawson] = lawson_iteration(r,z,f,grid,fvals,P(A),Q,sign_flag,nlawson,damping);
        if(maxerr_lawson<maxerr)
            % disp("Lawson decreased error by " + ((maxerr-maxerr_lawson)/maxerr) +"%");
            w = ww;
        end
    end

    % Compute the poles, residues, and zeroes of r(z)
    [pol, res, zer] = prz(z,f,w);
    % Repack r(z), enabling evaluation at support points
    numerator= @(x) prod(z.'-x,2) .* sum((f.*w).' ./ (x-z.'),2);
    denominator = @(x) prod(z.'-x,2) .* sum(w.' ./(x-z.'),2); 
    r = @(zz) reval(zz, z, f, w);
end



function [w,maxerr] = lawson_iteration(r,z,f,grid,fvals,A,Q,sign_flag,iteration_count,damping)
    wt = NaN; wtnew = ones(size(fvals));
    err = abs(r(grid)-fvals);
    for i = 1:iteration_count 
        % Update least-squares weights
        wt = wtnew; 
        wtnew = wt*(damping) + (1-damping)*wt.*err;
        % Solve for new weights
        w = new_weights(sqrt(wtnew).*A, @(M) M, Q, sign_flag);
        % Compute updated rational approximation
        numerator = @(x) sum((f.*w).' ./ (x-z.'),2);
        denominator = @(x) sum(w.' ./(x-z.'),2);
        r = @(x)  numerator(x) ./ denominator(x);
        %Compute errors
        err = abs(r(grid)-fvals);
    end
    maxerr = max(err);  
end

function [Q,P,d] = get_constraint(c,k,mm,nn,f,z,forZ3)
    d=0; Q = 0; a = 1;
    % Constrain numerator
    if(k > mm)
        d=k-mm-1;
        a = (f + forZ3*1).';
        [~,Q] = CoeffDivide(c,z,d,a); % TEMP solution
    % Constrain denominator
    elseif(k > nn) 
        d = k-nn-1;
        if(forZ3)
            a = (1-f).';
        end
        [~,Q] = CoeffDivide(c,z,d,a); % TEMP solution
    end
    P = @(A) A - (A*Q)*Q';
end

% MINIMIZES |Av|/|v| subject to v not in null(A);
function w = new_weights(A, P, Q, sign_flag)
    % Get right singular vectors not in null(A)
    [~,s,V] = svd(P(A),'econ','vector');
    % Throw away nearly-zero singular vectors
    keep = s>s(1)*1e-14;
    s = s(keep);
    V = V(:,keep);

    if(sign_flag)
        w = sum(V ./ (s.^2).', 2); % The 'sign improvement'
    else
        % w = V(:,end);
        mm = find( s == min(s) ); % Treat case of multiple min sing val
        nm = length(mm);
        w = V(:,mm)*ones(nm,1)/sqrt(nm); % Aim for non-sparse wt vector
    end
    % Project onto constraint
    w = w-Q*(Q'*w);
    w = w / norm(w);
end

% Updates a monic polynomial's coefficients given the next root.
% Uses the fact that p(z)(z-zk)=z*p(z)-zk*p(z).
% The first term shifts the coefficients, the second is a pointwise
% multiplication..
function c = update_polynomial_coeffs(c, zk, k, D)
    % Case 1: Polynomial is currently less than degree D
    if k < D+1
        c = -zk*c + [0 ; c(1:end-1)];
    % Case 2: Polynomial is at least degree D
    elseif k >= D+1
        c = -zk*[c(2:end);0] + c;
    end

end

% GREEDILY chooses next support point
function [zk,fk,j,maxerr] = select_next_support(grid, fvals, rk)
    % Compute the error over the grid
    err = abs(rk(grid)-fvals);
    % Greedy choice
    [maxerr,j] = max(err);

    % Get support point and function value there
    zk = grid(j);
    fk = fvals(j);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PRZ   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pol, res, zer] = prz(zj, fj, wj)
%   Compute poles, residues, and zeros of rational fun in barycentric form.

% Compute poles via generalized eigenvalue problem:
m = length(wj);
B = eye(m+1);
B(1,1) = 0;
E = [0 wj.'; ones(m, 1) diag(zj)];
pol = eig(E, B);
pol = pol(~isinf(pol));

% Compute residues via formula for res of quotient of analytic functions:
N = @(t) (1./(t-zj.')) * (fj.*wj);
Ddiff = @(t) -((1./(t-zj.')).^2) * wj;
res = N(pol)./Ddiff(pol);

% Compute zeros via generalized eigenvalue problem:
E = [0 (wj.*fj).'; ones(m, 1) diag(zj)];
zer = eig(E, B);
zer = zer(~isinf(zer));

end % End of PRZ.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   REVAL   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = reval(zz, zj, fj, wj)
%   Construct function handle to evaluate rational function in barycentric form.

zv = zz(:);                         % vectorize zz if necessary
CC = 1./(zv-zj.');                  % Cauchy matrix
r = (CC*(wj.*fj))./(CC*wj);         % vector of values

% Deal with input inf: r(inf) = lim r(zz) = sum(w.*f) / sum(w):
r(isinf(zv)) = sum(wj.*fj)./sum(wj);

% Deal with NaN:
ii = find(isnan(r));
for jj = 1:length(ii)
    if ( isnan(zv(ii(jj))) || ~any(zv(ii(jj)) == zj) )
        % r(NaN) = NaN is fine.
        % The second case may happen if r(zv(ii)) = 0/0 at some point.
    else
        % Clean up values NaN = inf/inf at support points.
        % Find the corresponding node and set entry to correct value:
        r(ii(jj)) = fj(zv(ii(jj)) == zj);
    end
end

% Reshape to input format:
r = reshape(r, size(zz));

end % End of REVAL.

