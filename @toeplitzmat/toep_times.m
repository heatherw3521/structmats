function h = toep_times(T, X)
%TOEP_TIMES Computes the Toeplitz-Matrix product TX efficiently
%   This is done through clever use of Fourier transform and
%   circulant matrices

    c1 = [T.tc ; 0; fliplr(T.tr(2:end).') ]';
    d = fft(c1)';
    % C = circulant(c1);
    
    R = [X;zeros(size(T,2), size(X,2))];
        
    y = ifft(d.*(fft(R)));
    h = y(1:size(T,1),:);
end

