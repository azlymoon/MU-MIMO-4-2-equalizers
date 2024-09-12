function [dw]=qam_demod_norm(rx,qam_seq)
    M=length(qam_seq);
    k=log2(M);
    dist_list=zeros(M,1);
    for jj=1:M
        dist_list(jj)=norm(rx - qam_seq(jj))^2;
    end
    pos = find(dist_list == min(dist_list)); 
    dw_sim=pos-1;
    dw=int2bit(dw_sim,k);
end