function evm_zf = calc_evm_zf(S, nOFDM, n, tx_qam, rx_fft_zf)
    evm_zf=zeros(S,nOFDM-2);

    for ss=1:S
        for jj=1:nOFDM-2
           evm_zf(ss,jj)=calculate_evm(tx_qam(ss,:,jj).',rx_fft_zf(:,jj,ss),n);
        end
    end
end
