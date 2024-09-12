function [dw_zf,dw_mmse] = demodulation(n, k, nOFDM, S, M, rx_fft_zf, rx_fft_mmse)
    dw_zf=zeros(n*k,nOFDM-2,S);
    dw_mmse=zeros(n*k,nOFDM-2,S);
    for ss=1:S
        for jj=1:nOFDM-2
            for ii=0:n-1
                %Choose rx_fft_zf or rx_fft_mmse
                %dw(ii*k+1:ii*k+k,jj)=qam_demod(rx_fft_zf(ii+1,jj),qam_seq_gray,Es);
                %dw(ii*k+1:ii*k+k,jj)=qam_demod_norm(rx_fft_zf(ii+1,jj),qam_seq_gray);
                dw_zf(ii*k+1:ii*k+k,jj,ss)=qamdemod(rx_fft_zf(ii+1,jj,ss),M,'gray',OutputType='bit',UnitAveragePower=true);
                dw_mmse(ii*k+1:ii*k+k,jj,ss)=qamdemod(rx_fft_mmse(ii+1,jj,ss),M,'gray',OutputType='bit',UnitAveragePower=true);
            end
        end
    end
end

