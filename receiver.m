function rx_fft = receiver(rx_awgn, n,nOFDM,nFFT,padding_len,CP_len,delay, L)
    rx_fft = zeros(n,nOFDM,L);
    for l=1:L
        rx_fft(:,:,l) = OFDM_demod(rx_awgn(l,:),n,nOFDM,nFFT,padding_len,CP_len,delay);
    end
end

