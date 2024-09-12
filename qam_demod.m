function [dw]=qam_demod(rx,qam_seq,Es)

    y=sqrt(Es)*rx;

    M=length(qam_seq); 

    rIdx = 2*floor(real(y)/2)+1;
    if abs(rIdx)>(sqrt(M)-1)
        if real(y)>0
            rIdx=sqrt(M)-1;
        else
            rIdx=-sqrt(M)+1;
        end
    end

    iIdx = 2*floor(imag(y)/2)+1;

    if abs(iIdx)>(sqrt(M)-1)
        if real(y)>0
            iIdx=sqrt(M)-1;
        else
            iIdx=-sqrt(M)+1;
        end
    end

    dw=complex(rIdx,iIdx);
    pos=find(qam_seq==dw);
    dw=int2bit(pos-1,log2(M));
end