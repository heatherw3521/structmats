function h = ptwise_apply(f, V1, V2)
%PTWISE_APPLY Pointwise application of f to (V1,V2), where at least on of
%   V1,V2 is a VANDERMAT. Unfortunately, I don't think there is much
%   structure here to take advantage of. The only pointwise operation that
%   retains vandermonde structure is V .^ a, where a is a scalar. This is
%   handled in a separate file

end

