function [t, ts, out] = mp_timeit(f, reps, warmup)
%MP_TIMEIT  Median wall time of f() over `reps` runs after `warmup` runs.
%   (Octave has no timeit; tic/toc is used for both.)  out = last output.
if nargin < 2, reps = 3; end
if nargin < 3, warmup = 1; end
for k = 1:warmup, out = f(); end %#ok<NASGU>
ts = zeros(reps,1);
for k = 1:reps
    t0 = tic; out = f(); ts(k) = toc(t0);
end
t = median(ts);
end
