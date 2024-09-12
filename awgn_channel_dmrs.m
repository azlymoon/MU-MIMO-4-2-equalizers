function rx=awgn_channel_dmrs(tx_ifft_seq,nFFT,CP_len,sym_1,sym_2,noise_var,delay)
    %%%%%AWGN channel%%%%%%
    noise = (randn(size(tx_ifft_seq)) + 1i * randn(size(tx_ifft_seq))) / sqrt(2);
    noise((nFFT + CP_len) * (sym_1 - 1) + 1 + delay:(nFFT + CP_len) * sym_1 + delay) = complex(zeros(1, nFFT + CP_len));
    noise((nFFT + CP_len) * (sym_2 - 1) + 1 + delay:(nFFT + CP_len) * sym_2 + delay) = complex(zeros(1, nFFT + CP_len));
    rx = tx_ifft_seq + sqrt(noise_var) * noise;
end