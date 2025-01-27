clear all
fe = 10000;       %the number of collected statistics at each noise var value

%Modulation parameters
M = 16;      % Modulation order
k = log2(M); % Number of bits per symbol
[qam_seq_bin,qam_list] = get_sequence(M);

%Average power of m-ary symbol for QAM
Es=0;
for ii = 1:M
    Es=Es+(qam_list(ii, 1) + qam_list(ii, 2))^2;
end
Es=Es/M;

%OFDM parameters
nOFDM = 14;                % Number of OFDM symbols
nFFT = 2048;             % Number of FFT bins
rb = 100;                % Number of RB (resource blocks)
n = rb*12;               % Number of QAM symbols in OFDM symbol
padding_len = (nFFT - n) / 2;  % length of padding zeros elements before ifft (half of them)
CP_len = 0.125 * nFFT;       % length of cyclic prefix (guard band)

%Fading channel parameters
fs = 30.72e6;            % Sample rate in Hz from LTE

%Number of antennas
L = 8;

%Number of streams
S = 4;

delay=7;

%LTE model
model = struct( ...
    DelayProfile = "EVA", ...
    NRxAnts=1, ...
    DopplerFreq = 5, ...
    MIMOCorrelation = "Low", ...
    Seed = 0, ...
    ModelType = "Dent", ...
    NTerms = 16, ...
    NormalizeTxAnts = "On", ...
    NormalizePathGains = "On" ...
    );
model.SamplingRate = fs;
model.InitTime = 0;

%DM-RS sequnces
ue.NCellID = 3;
ue.NSubframe = 0;

chs.PRBSet = (0:49).';

dmrs_1 = ltePUSCHDRS(ue,chs);

ue.NCellID = 10;
dmrs_2 = ltePUSCHDRS(ue,chs);

sym_1 = 4; sym_2 = 11;

%Parameters for simulations
EbNoArray = -2:2:8;
Bit_Err_zf = zeros(1, length(EbNoArray));
Bit_Err_mmse = zeros(1, length(EbNoArray));
Bit_Err_Streams_zf = zeros(S, length(EbNoArray));
Bit_Err_Streams_mmse = zeros(S, length(EbNoArray));
BER_theor = zeros(1, length(EbNoArray));
Sim_Err = zeros(1, length(EbNoArray));
SER_theor = zeros(1, length(EbNoArray));
Itr = zeros(1, length(EbNoArray));
evm_zf_avg = zeros(1, length(EbNoArray));
evm_mmse_avg = zeros(1, length(EbNoArray));
evm_terminal_zf_avg = zeros(S, length(EbNoArray));
evm_terminal_mmse_avg = zeros(S, length(EbNoArray));

%Converting EsNO from EbNo
EsNoArray = EbNoArray + 10*log10(k) + 10*log10(n/nFFT);

count = 0;
for i = 1:length(EbNoArray)
    while count < fe
        EsNo = EsNoArray(i);
        noise_var = 1/(exp(log(10) * EsNo/10)); %QAM with unit average power

        %%%%%Transmitter%%%%%%

        %Generating data subframe w/o dm-rs positions
        data = randi([0 1], n*k, nOFDM-2, S);

        [tx_ifft_seq_L, tx_qam] = transmitter(data, n, nOFDM, S, nFFT, CP_len, delay, M,sym_1,sym_2,dmrs_1,dmrs_2,padding_len, L);


        %%%%% Сhannel %%%%%%
        rx_awgn = channel(tx_ifft_seq_L, nOFDM, S, nFFT, CP_len, delay, sym_1,sym_2, L, model, noise_var);


        %%%%%Receiver%%%%%%
        rx_fft = receiver(rx_awgn, n,nOFDM,nFFT,padding_len,CP_len,delay, L);


        %Equalizer
        %Channel response
        h_channel = channel_response(rx_fft, sym_1, sym_2, L, n, S, dmrs_1, dmrs_2);


        %ZF equalizer
        rx_fft_zf = zf_equalizer(h_channel, rx_fft, n, S, L, nOFDM, sym_1, sym_2);


        %MMSE equalizer
        rx_fft_mmse = mmse_equalizer(h_channel, rx_fft, n, S, L, nOFDM, sym_1, sym_2, noise_var);


        %EVM
        evm_zf = calc_evm_zf(S, nOFDM, n, tx_qam, rx_fft_zf);
        evm_mmse = calc_evm_mmse(S, nOFDM, n, tx_qam, rx_fft_mmse);

        evm_zf_avg(i) = evm_zf_avg(i) + mean(10*log10(evm_zf),'all');
        evm_mmse_avg(i) = evm_mmse_avg(i) + mean(10*log10(evm_mmse),'all');

        for ss=1:S
            evm_terminal_zf_avg(ss, i) = evm_terminal_zf_avg(i)+mean(10*log10(evm_zf(ss, :)),'all');
            evm_terminal_mmse_avg(ss, i) = evm_terminal_mmse_avg(i)+mean(10*log10(evm_mmse(ss, :)),'all');
        end


        %Demodulation
        [dw_zf, dw_mmse] = demodulation(n, k, nOFDM, S, M, rx_fft_zf, rx_fft_mmse);

        %Count SER and BER
        bit_err_zf=0;
        bit_err_mmse=0;

        for ss = 1:S
            for jj = 1:nOFDM-2
                for ii = 0:n-1
                    bit_err_cur_zf = nnz(data(ii*k+1 : ii*k+k, jj, ss) - dw_zf(ii*k+1 : ii*k+k, jj, ss));
                    bit_err_cur_mmse = nnz(data(ii*k+1 : ii*k+k, jj, ss) - dw_mmse(ii*k+1 : ii*k+k, jj, ss));
                    if bit_err_cur_zf > 0
                        count = count + 1;
                        Sim_Err(i) = Sim_Err(i) + 1;
                    end
                    bit_err_zf = bit_err_zf + bit_err_cur_zf;
                    bit_err_mmse = bit_err_mmse + bit_err_cur_mmse;
                end
            end
            Bit_Err_Streams_zf(ss, i) = Bit_Err_Streams_zf(ss, i) + bit_err_zf;
            Bit_Err_Streams_mmse(ss, i) = Bit_Err_Streams_mmse(ss, i) + bit_err_mmse;
        end
        Bit_Err_zf(i) = Bit_Err_zf(i) + bit_err_zf;
        Bit_Err_mmse(i) = Bit_Err_mmse(i) + bit_err_mmse;
        Itr(i) = Itr(i) + 1;
    end
    count = 0;
    %BER_theor(i)=(1/k)*2*(sqrt(M)-1)/sqrt(M)*erfc(sqrt(k*0.1*(10.^(EbNoArray(i)/10))));
    %SER_theor(i)=2*(sqrt(M)-1)/sqrt(M)*erfc(sqrt(0.1*(10.^(EsNoArray(i)/10))));
    [BER_theor(i), SER_theor(i)]=berfading(EbNoArray(i),'qam',M,S,L);

    disp("-------------");
    disp(Bit_Err_mmse./Itr./(S*n*k*(nOFDM-2)));
    disp(Bit_Err_zf./Itr./(S*n*k*(nOFDM-2)));
end

BER_zf = Bit_Err_zf./Itr./(S*n*k*(nOFDM-2));
BER_mmse = Bit_Err_mmse./Itr./(S*n*k*(nOFDM-2));
BER_Terminal_zf = Bit_Err_Streams_zf./Itr./(n*k*(nOFDM-2));
BER_Terminal_mmse = Bit_Err_Streams_mmse./Itr./(n*k*(nOFDM-2));
SER = Sim_Err./Itr./(S*n*(nOFDM-2));
EVM_zf = evm_zf_avg./Itr;
EVM_mmse = evm_mmse_avg./Itr;
EVM_Terminal_zf = evm_terminal_zf_avg./Itr;
EVM_Terminal_mmse = evm_terminal_mmse_avg./Itr;
