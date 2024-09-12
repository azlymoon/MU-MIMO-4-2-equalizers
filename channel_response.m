function h_channel = channel_response(rx_fft, sym_1, sym_2, L, n, S, dmrs_1, dmrs_2)
    h_channel = zeros([L, n, S]);
    
    for l = 1:L
        h_channel(l, :, 1) = rx_fft(:,sym_1,l)./dmrs_1;
        h_channel(l, :, 2) = rx_fft(:,sym_2,l)./dmrs_2;
    end
end
