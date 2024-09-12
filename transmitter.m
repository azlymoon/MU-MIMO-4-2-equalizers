function [tx_ifft_seq_L, tx_qam] = transmitter(data, n, nOFDM, S, nFFT, CP_len, delay, M,sym_1,sym_2,dmrs_1,dmrs_2,padding_len, L)
        tx_ifft_seq = zeros(S,nOFDM*(nFFT+CP_len)+delay);
        tx_qam = zeros(S,n,(nOFDM-2));
        for s=1:S
            [tx_ifft_seq(s,:),tx_qam(s,:,:)]=OFDM_modulation(data(:,:,s),n,M,nOFDM,sym_1,sym_2,dmrs_1,dmrs_2,nFFT,padding_len,CP_len,delay,s);
        end

        % Duplicate the IFFT sequence for L antennas
        % tx_ifft_seq_L has L*S rows
        tx_ifft_seq_L = repelem(tx_ifft_seq, L, 1);
end

