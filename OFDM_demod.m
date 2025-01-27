function rx_fft = OFDM_demod(rx, n, nOFDM, nFFT, padding_len, CP_len, delay)
        norm = nFFT / sqrt(n);
        %Serial to parallel
        rx_cp = reshape(rx(delay+1:end), [nFFT+CP_len, nOFDM]);

        %CP removal
        rx = rx_cp(CP_len+1:end, :);

        %FFT step
        rx_fft_padding = fft(rx, nFFT, 1)./norm;
        rx_fft_padding = fftshift(rx_fft_padding, 1);

        rx_fft = rx_fft_padding(padding_len+1:(end-padding_len), :);
end