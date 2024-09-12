function [qam_seq_bin, qam_list] = get_sequence(M)

    qam_seq=[sqrt(M)-1:-2:1 -1:-2:-(sqrt(M)-1)];
    qam_list=table2array(combinations(flip(qam_seq),qam_seq));

    qam_complex=zeros(M,1);
    for ii=1:length(qam_list)
        qam_complex(ii)=qam_list(ii,1)+i*qam_list(ii,2);
    end

    qam_seq_bin=qam_complex;
end