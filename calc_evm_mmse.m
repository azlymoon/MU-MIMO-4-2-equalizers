function evm_mmse = calc_evm_mmse(S, nOFDM, n, tx_qam, rx_fft_mmse)
    evm_mmse=zeros(S,nOFDM-2);
    
    for ss=1:S
        for jj=1:nOFDM-2
            evm_mmse(ss,jj)=calculate_evm(tx_qam(ss,:,jj).',rx_fft_mmse(:,jj,ss),n);
        end
    end
end

