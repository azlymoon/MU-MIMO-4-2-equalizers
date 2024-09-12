% Функция для вычисления EVM
function evm = calculate_evm(tx_qam, rx_fft_eq, n)
    % Вычисление EVM для ZF Equalizer
    evm = (1 / (n)) * sum((imag(tx_qam) - imag(rx_fft_eq)).^2 + (real(tx_qam) - real(rx_fft_eq)).^2, 'all');
end