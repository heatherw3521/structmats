function [t_first, t_repeat, x] = mp_time_solve(H, b, reps)
%MP_TIME_SOLVE  Wall times of H\b without and with stored factors.
%   H keeps its factors after a solve (@hss/hss.m), so a plain repeated
%   timing of H\b would only measure the stored-factor solve. This times both:
%     t_first   median over reps of  clearfactors(H); H\b   (factor + solve)
%     t_repeat  median over reps of  H\b  with the factors in place
%   One untimed warm-up solve comes first.
if nargin < 3, reps = 3; end
clearfactors(H); x = H \ b;                       % warm-up
tf = zeros(reps, 1);
for k = 1:reps
    clearfactors(H);
    t0 = tic; x = H \ b; tf(k) = toc(t0);
end
tr = zeros(reps, 1);
for k = 1:reps
    t0 = tic; x = H \ b; tr(k) = toc(t0);
end
t_first = median(tf);
t_repeat = median(tr);
end
