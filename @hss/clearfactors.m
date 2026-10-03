function H = clearfactors(H)
%CLEARFACTORS  Free the factorizations stored with H.
%   H\b, minnorm and tikhonov keep their factors with H (about the memory of H
%   itself per stored factorization) so that later calls only solve.
%   clearfactors(H) drops them; the next solve factors again. Copies of H that
%   share the factors (H2 = H, unmodified) lose them too.
if ~isempty(H.factorcache)
    H.factorcache.reset();
end
end
