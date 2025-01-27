function rx_fft_mmse = mmse_equalizer(h_channel, rx_fft, n, S, L, nOFDM, sym_1, sym_2, noise_var)
    W = zeros([n, S, L]);

    for ii = 1:n
       W(ii, :, :) = ((squeeze(h_channel(:, ii, :))' * squeeze(h_channel(:, ii, :)) + noise_var))\squeeze(h_channel(:, ii, :))';
    end

    % W(n/2+1:end, :, :) = W(1:n/2, :, :);
    
    rx_fft_eq_1 = zeros(n,nOFDM,S);
     
    for jj = 1:14
        for ii = 1:n
            rx_fft_eq_1(ii,jj,:) = squeeze(W(ii, :, :)) * squeeze(rx_fft(ii, jj, :));
        end
    end

    rx_fft_mmse = rx_fft_eq_1(:, [1:sym_1-1, sym_1+1:sym_2-1, sym_2+1:end], :);
end

