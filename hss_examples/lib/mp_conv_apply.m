function y = mp_conv_apply(K, x)
%MP_CONV_APPLY  Exact y = A x for a mp_kernel_conv operator (via conv), any size.
c = conv(x(:), K.g);                      % c(k) = sum_j x_j g(k - j - w) (1-based shift)
r = K.stride*(0:K.m-1)' + K.offset;
y = c(r + K.w + 1);
end
