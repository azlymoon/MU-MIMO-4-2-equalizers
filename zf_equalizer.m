function rx_fft_zf = zf_equalizer(h_channel, rx_fft, n, S, L, nOFDM, sym_1, sym_2)
    W = zeros([n, S, L]);
            
    for ii = 1:n
        W(ii, :, :) = pinv(squeeze(h_channel(:, ii, :)));
    end

    rx_fft_eq=zeros(n,nOFDM,S);
    
    for jj = 1:14
        for ii=1:n
            rx_fft_eq(ii,jj,:)=squeeze(W(ii,:,:))*squeeze(rx_fft(ii,jj,:));
        end
    end

    rx_fft_zf=rx_fft_eq(:,[1:sym_1-1,sym_1+1:sym_2-1,sym_2+1:end],:);
end

