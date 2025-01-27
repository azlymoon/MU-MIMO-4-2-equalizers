function rx_awgn = channel(tx_ifft_seq_L, nOFDM, S, nFFT, CP_len, delay, sym_1,sym_2, L, model, noise_var)
    %Fading
    rx_fading_lte = zeros(size(tx_ifft_seq_L));
    for l = 1:L*S
        rx_fading_lte(l,:) = lteFadingChannel(model, reshape(tx_ifft_seq_L(l,:), [nOFDM*(nFFT+CP_len) + delay, 1]));
    end
    
    %Superposition on L antennas
    rx_superposition = zeros(L, size(rx_fading_lte, 2));
    for s = 1:S
        rx_superposition = rx_superposition + rx_fading_lte(1+L*(s-1) : s*L, : );
    end
    
    rx_awgn=zeros(L, size(rx_superposition,2));
    for l = 1:L
        rx_awgn(l, :) = awgn_channel_dmrs(rx_superposition(l, :), nFFT, CP_len, sym_1, sym_2, noise_var, delay);
    end
end
