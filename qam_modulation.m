function [tx,qam_seq_gray] = qam_modulation(M,data,qam_seq,Es)
    symbols = bit2int(data,length(data));
    [~, seq_gray] = comm.internal.qam.getGrayMapping(M, symbols);
    qam_seq_gray = qam_seq(seq_gray+1);
    tx = qam_seq_gray(symbols+1)/sqrt(Es);
end