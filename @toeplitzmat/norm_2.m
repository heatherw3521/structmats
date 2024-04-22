function h = norm_2(T)
%NORM_2 Efficiently computes the 2-norm for Toeplitz matrices.
%   The 2-norm is equal to the largest singular value.
%
%   We solve this problem by using the power method, taking advantage of
%   fast toeplitz-vector products.

    m = size(T,1);
    n = size(T,2);

    if(n < m) % Let's deal with smaller vectors
        h = norm_2(T.');
    else % We use T.'*T
        d = 1; lam = 1;
        Tstar = T';
        v = rand(n,1);
        while 1
            v = Tstar * (T*v); v = v / norm(v);
            lam_new = v.' * (T.'*(T*v));
            d = abs(lam - lam_new);
            lam = lam_new;
            if (d/lam < 1e-12)
                break;
            end
        end
        h = sqrt(lam);
    end
end

