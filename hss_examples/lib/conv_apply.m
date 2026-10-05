function y = conv_apply(K, x)
%CONV_APPLY  Exact y = A x for a kernel_conv operator (via conv), any size.
c = conv(x(:), K.g);                      % c(k) = sum_j x_j g(k - j - w) (1-based shift)
r = K.stride*(0:K.m-1)' + K.offset;
y = c(r + K.w + 1);
end
