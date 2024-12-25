function h_channel = channel_response(rx_fft, sym_1, sym_2, L, n, S, dmrs_1, dmrs_2)
    h_channel = zeros([L, n, S]);
    
    if S == 2
        for l = 1:L
            h_channel(l, :, 1) = rx_fft(:,sym_1,l)./dmrs_1;
            h_channel(l, :, 2) = rx_fft(:,sym_2,l)./dmrs_2;
        end
    elseif S == 4
        for l = 1:L
            even_indices = (2:2:n)';
            odd_indices = (1:2:n)';

            h_channel(l, even_indices, 1) = rx_fft(even_indices, sym_1, l) ./ dmrs_1(even_indices);
            h_channel(l, odd_indices, 1) = interp1(even_indices, h_channel(l, even_indices, 1), odd_indices, 'linear', 'extrap');
                
            h_channel(l, odd_indices, 2) = rx_fft(odd_indices, sym_1, l) ./ dmrs_1(odd_indices);
            h_channel(l, even_indices, 2) = interp1(odd_indices, h_channel(l, odd_indices, 2), even_indices, 'linear', 'extrap');

            h_channel(l, even_indices, 3) = rx_fft(even_indices, sym_2, l) ./ dmrs_2(even_indices);
            h_channel(l, odd_indices, 3) = interp1(even_indices, h_channel(l, even_indices, 3), odd_indices, 'linear', 'extrap');

            h_channel(l, odd_indices, 4) = rx_fft(odd_indices, sym_2, l) ./ dmrs_2(odd_indices);
            h_channel(l, even_indices, 4) = interp1(odd_indices, h_channel(l, odd_indices, 4), even_indices, 'linear', 'extrap');
        end
    end
