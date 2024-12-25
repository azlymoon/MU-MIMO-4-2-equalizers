function [tx_ifft_seq, tx_qam] = OFDM_modulation(data, n, M, nOFDM, sym_1, sym_2, dmrs_1, dmrs_2, nFFT, padding_len, CP_len, delay, S, current_s)
    k = log2(M);
    
    % Модуляция
    tx_qam = zeros(n, nOFDM-2);
    for jj = 1:nOFDM-2
        for ii = 0:n-1
            tx_qam(ii+1,jj) = qammod(data(ii*k+1:ii*k+k,jj), M, 'gray', 'InputType', 'bit', 'UnitAveragePower', true);
        end
    end

    % Добавление dm-rs сигнала на позиции sym_1 и sym_2
    tx_qam_dmrs = zeros(n, nOFDM);

    % Заполнение tx_qam_dmrs данными
    tx_qam_dmrs(:, [1:sym_1-1, sym_1+1:sym_2-1, sym_2+1:end]) = tx_qam;

    if S == 1
        tx_qam_dmrs(:, sym_1) = dmrs_1;
        tx_qam_dmrs(:, sym_2) = dmrs_2;

    elseif S == 2
        if current_s == 1
            tx_qam_dmrs(:, sym_1) = dmrs_1;
            tx_qam_dmrs(:, sym_2) = 0;
        elseif current_s == 2
            tx_qam_dmrs(:, sym_1) = 0;
            tx_qam_dmrs(:, sym_2) = dmrs_2;
        end

    elseif S == 4
        even_indices = (2:2:n)';
        odd_indices = (1:2:n)';
        column_sym_1 = zeros(size(tx_qam_dmrs(:, sym_1)));
        column_sym_2 = zeros(size(tx_qam_dmrs(:, sym_2)));
        if current_s == 1
            column_sym_1(even_indices) = dmrs_1(even_indices);
        elseif current_s == 2
            column_sym_1(odd_indices) = dmrs_1(odd_indices);
        elseif current_s == 3
           column_sym_2(even_indices) = dmrs_2(even_indices);
        elseif current_s == 4
            column_sym_2(odd_indices) = dmrs_2(odd_indices);
        end

        tx_qam_dmrs(:, sym_1) = column_sym_1;
        tx_qam_dmrs(:, sym_2) = column_sym_2;
    end

    % iFFT шаг с нулевым заполнением
    tx_padding = [complex(zeros(padding_len, nOFDM)); tx_qam_dmrs; complex(zeros(padding_len, nOFDM))];
    tx_padding = ifftshift(tx_padding, 1);
    norm = nFFT / sqrt(n);
    tx_ifft = norm * ifft(tx_padding, nFFT, 1);

    % Добавление CP
    tx_ifft_cp = [tx_ifft(end-CP_len+1:end, :); tx_ifft];

    % Параллельный в последовательный
    tx_ifft_seq = reshape(tx_ifft_cp, [1, nOFDM * (nFFT + CP_len)]);
    tx_ifft_seq = [tx_ifft_seq, zeros(1, delay)];
end
