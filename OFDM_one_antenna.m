clear all
fe=20000;       %the number of collected statistics at each noise var value

%Modulation parameters
M = 4;      % Modulation order
k = log2(M); % Number of bits per symbol
[qam_seq_bin,qam_list]=get_sequence(M);

%Average power of m-ary symbol for QAM
Es=0;
for ii = 1:M
    Es=Es+(qam_list(ii,1)+qam_list(ii,2))^2;
end
Es=Es/M;

%OFDM parameters
nOFDM=14;                % Number of OFDM symbols   
nFFT = 2048;             % Number of FFT bins
rb = 10;                % Number of RB (resource blocks)
n = rb*12;               % Number of QAM symbols in OFDM symbol
padding_len=(nFFT-n)/2;  % length of padding zeros elements before ifft (half of them)
CP_len=0.125*nFFT;       % length of cyclic prefix (guard band)

%Fading channel parameters
fs = 30.72e6;            % Sample rate in Hz from LTE

%Number of antennas
L = 4;

%LTE model
model = struct(DelayProfile="EVA",NRxAnts=1, ...
    DopplerFreq = 5,MIMOCorrelation="Low", ...
    Seed=0,ModelType="Dent", ...
    NTerms=16,NormalizeTxAnts="On", ...
    NormalizePathGains="On");
model.SamplingRate=fs;
model.InitTime= 0;

%DM-RS sequnces
ue.NCellID = 3;
ue.NSubframe = 0;

chs.PRBSet = (0:4).';
%%chs.PRBSet = (0:49).';

dmrs_1 = ltePUSCHDRS(ue,chs);

ue.NCellID = 10;
dmrs_2 = ltePUSCHDRS(ue,chs);

sym_1=4; sym_2=11;

%Parameters for simulations
EbNoArray=0:2:12;
Bit_Err_zf=zeros(1,length(EbNoArray));
Bit_Err_mmse=zeros(1,length(EbNoArray));
BER_theor=zeros(1,length(EbNoArray));
Sim_Err=zeros(1,length(EbNoArray));
SER_theor=zeros(1,length(EbNoArray));
Itr=zeros(1,length(EbNoArray));
evm_zf_avg=zeros(1,length(EbNoArray));
evm_mmse_avg=zeros(1,length(EbNoArray));

%Converting EsNO from EbNo
EsNoArray=EbNoArray+10*log10(k);%+10*log10(n/nFFT);


% IRC: Разложение Холецкого ковариационной матрицы
K_ni = [1 0.9 0.9 0.81; 0.9 1 0.81 0.9; 0.9 0.81 1 0.9; 0.81 0.9 0.9 1];
L_matrix = chol(K_ni, 'lower');


count=0;
for i=1:length(EbNoArray)
    while count<fe
        EsNo=EsNoArray(i);
        noise_var = 1/(exp(log(10) * EsNo/10)); %QAM with unit average power

        %%%%%Transmitter%%%%%%

        %Generating data subframe w/o dm-rs positions
        data = randi([0 1],n*k,nOFDM-2);

        %Modulation
        tx_qam=zeros(n,nOFDM-2);
        for jj=1:nOFDM-2
            for ii=0:n-1
                tx_qam(ii+1,jj)=qammod(data(ii*k+1:ii*k+k,jj),M,'gray',InputType='bit',UnitAveragePower=true);
            end
        end

        %Add dm-rs signal on positions sym_1 and sym_2
        tx_qam_dmrs=zeros(n,nOFDM);

        tx_qam_dmrs(:,[1:sym_1-1,sym_1+1:sym_2-1,sym_2+1:end])=tx_qam;
        tx_qam_dmrs(:,sym_1)=dmrs_1;
        tx_qam_dmrs(:,sym_2)=dmrs_2;

        %iFFt step with zero padding
        tx_padding=[complex(zeros(padding_len,nOFDM)); tx_qam_dmrs; complex(zeros(padding_len,nOFDM))];
        tx_padding=ifftshift(tx_padding,1);
        norm=sqrt(nFFT);%/sqrt(n);
        tx_ifft=norm*ifft(tx_padding,nFFT,1);

        %adding CP
        tx_ifft_cp=[tx_ifft(end-CP_len+1:end,:); tx_ifft];

        delay=7;
        %Parallel to serial
        tx_ifft_seq=reshape(tx_ifft_cp,[1,nOFDM*(nFFT+CP_len)]);
        tx_ifft_seq=[tx_ifft_seq, zeros(1,delay)];

        % Duplicate the IFFT sequence for L antennas
        tx_ifft_seq_L = repmat(tx_ifft_seq, L, 1);

        %%%%%Fading channel%%%%%%
        rx_fading_lte = zeros(length(tx_ifft_seq),L);
        for l = 1:L
           [rx_fading_lte(:, l), info] = lteFadingChannel(model, reshape(tx_ifft_seq_L(l, :), [numel(tx_ifft_seq), 1]));
        end

        %rx_fading_lte = tx_ifft_seq_L.';
        %%%%%Shift offset%%%%%%%%
        %imp_res=[zeros(1,delay), 1, zeros(1,200)];
        %rx_fading_phase=conv(tx_ifft_seq,imp_res);
        %rx_fading_phase=reshape(rx_fading_phase(1:length(tx_ifft_seq)),[1,length(tx_ifft_seq)]);

        %%%%%AWGN channel%%%%%%
        noise = zeros(size(rx_fading_lte));
        % rx = zeros(size(rx_fading_lte));
        for l = 1:L
            noise(:, l) = (randn(size(rx_fading_lte(:, l))) + 1i * randn(size(rx_fading_lte(:, l)))) / sqrt(2);
            noise((nFFT + CP_len) * (sym_1 - 1) + 1 + delay:(nFFT + CP_len) * sym_1 + delay, l) = complex(zeros(1, nFFT + CP_len));
            noise((nFFT + CP_len) * (sym_2 - 1) + 1 + delay:(nFFT + CP_len) * sym_2 + delay, l) = complex(zeros(1, nFFT + CP_len));
        end
        % Добавляем L_matrix в AWGN
        rx = rx_fading_lte + sqrt(noise_var) * noise * L_matrix;

        %%%%%Receiver%%%%%%
        %Seral to parallel
        rx_cp = zeros([size(tx_ifft_cp), L]);
        for l = 1:L
            rx_cp(:, :, l)=reshape(rx(delay+1:end,l),[nFFT+CP_len,nOFDM]);
        end

        %CP removal
        rx = zeros([size(tx_ifft), L]);
        for l = 1:L
            rx(:, :, l)=rx_cp(CP_len+1:end,:, l);
        end

        %FFT step
        rx_fft_padding = zeros([size(tx_ifft), L]);
        for l = 1:L
            rx_fft_padding(:, :, l)=fft(rx(:,:,l),nFFT,1)./norm;
            rx_fft_padding(:, :, l)=fftshift(rx_fft_padding(:, :, l),1);
        end
        
        rx_fft = zeros([size(tx_qam_dmrs), L]);
        for l = 1:L
            rx_fft(:, :, l)=rx_fft_padding(padding_len+1:(end-padding_len),:, l);
        end

        %Equalizer
        %Channel response
        h_channel = zeros([L, n, 2]);
        for l = 1:L
            h_channel(l, :, 1) = rx_fft(:,sym_1,l)./dmrs_1;
            h_channel(l, :, 2) = rx_fft(:,sym_2,l)./dmrs_2;
        end

        %ZF equalizer
        W = zeros([n, L, 2]);
        for jj = 1:2
            for ii = 1:n
                W(ii, :, jj) = pinv(h_channel(:, ii, jj));
            end
        end

        rx_fft_eq=zeros(n,nOFDM);
        for j=1:7
            for ii=1:n
                rx_fft_eq(ii,j)=W(ii,:,1)*squeeze(rx_fft(ii,j,:));
                rx_fft_eq(ii,j+7)=W(ii,:,2)*squeeze(rx_fft(ii,j+7,:));
            end
        end    

        % MMSE equalizer with IRC
        % ------------
        % IRC: Применение фильтра предварительного подавления шум
        h_channel_filtered = zeros(size(h_channel));
        rx_fft_filtered = zeros(size(rx_fft));
        W = zeros([n, L, 2]);

        for dmrs_symbol = 1:2
            for ii = 1:n
                h_channel_filtered(:, ii, dmrs_symbol) = L_matrix \ h_channel(:, ii, dmrs_symbol);
            end
        end

        for ii = 1:n
            for jj = 1:nOFDM
                rx_fft_filtered(ii, jj, :) = L_matrix \ squeeze(rx_fft(ii, jj, :));
            end
        end

        for jj = 1:2
            for ii = 1:n
                W(ii, :, jj) = ((h_channel_filtered(:, ii, jj)' * h_channel_filtered(:, ii, jj) + noise_var))\h_channel_filtered(:, ii, jj)';
            end
        end

        rx_fft_eq_1=zeros(n,nOFDM);

        for j=1:7
            for ii=1:n
                rx_fft_eq_1(ii,j)=W(ii,:,1)*squeeze(rx_fft_filtered(ii,j,:));
                rx_fft_eq_1(ii,j+7)=W(ii,:,2)*squeeze(rx_fft_filtered(ii,j+7,:));
            end
        end
        % ------------

        % MMSE equalizer without IRC
        % ------------
        % W = zeros([n, L, 2]);
        % for jj = 1:2
        %     for ii = 1:n
        %         W(ii, :, jj) = ((h_channel(:, ii, jj)' * h_channel(:, ii, jj) + noise_var))\h_channel(:, ii, jj)';
        %     end
        % end
        % 
        % rx_fft_eq_1=zeros(n,nOFDM);
        % 
        % for j=1:7
        %     for ii=1:n
        %         rx_fft_eq_1(ii,j)=W(ii,:,1)*squeeze(rx_fft(ii,j,:));
        %         rx_fft_eq_1(ii,j+7)=W(ii,:,2)*squeeze(rx_fft(ii,j+7,:));
        %     end
        % end
        % ------------

        rx_fft_zf=rx_fft_eq(:,[1:sym_1-1,sym_1+1:sym_2-1,sym_2+1:end]);
        rx_fft_mmse=rx_fft_eq_1(:,[1:sym_1-1,sym_1+1:sym_2-1,sym_2+1:end]);

        evm_zf=zeros(1,nOFDM-2);
        evm_mmse=zeros(1,nOFDM-2);
        for jj=1:nOFDM-2
            evm_zf(jj)=calculate_evm(tx_qam(:,jj),rx_fft_zf(:,jj),n);
            evm_mmse(jj)=calculate_evm(tx_qam(:,jj),rx_fft_mmse(:,jj),n);
        end

        evm_zf_avg(i)=evm_zf_avg(i)+mean(10*log10(evm_zf));
        evm_mmse_avg(i)=evm_mmse_avg(i)+mean(10*log10(evm_mmse));

        dw_zf=zeros(n*k,nOFDM-2);
        dw_mmse=zeros(n*k,nOFDM-2);
        for jj=1:nOFDM-2
            for ii=0:n-1
                %Choose rx_fft_zf or rx_fft_mmse
                dw_mmse(ii*k+1:ii*k+k,jj)=qamdemod(rx_fft_mmse(ii+1,jj),M,'gray',OutputType='bit',UnitAveragePower=true);
                dw_zf(ii*k+1:ii*k+k,jj)=qamdemod(rx_fft_zf(ii+1,jj),M,'gray',OutputType='bit',UnitAveragePower=true);
            end
        end

        %Count SER and BER
        bit_err_zf=0;
        bit_err_mmse=0;
        for jj=1:nOFDM-2 
            for ii=0:n-1
                bit_err_cur_zf = nnz(data(ii*k+1:ii*k+k,jj)-dw_zf(ii*k+1:ii*k+k,jj));
                bit_err_cur_mmse = nnz(data(ii*k+1:ii*k+k,jj)-dw_mmse(ii*k+1:ii*k+k,jj));
                if bit_err_cur_mmse>0
                    count=count+1;
                    Sim_Err(i)=Sim_Err(i)+1;
                end
                bit_err_zf=bit_err_zf+bit_err_cur_zf;
                bit_err_mmse=bit_err_mmse+bit_err_cur_mmse;
            end
        end
        Bit_Err_zf(i)=Bit_Err_zf(i)+bit_err_zf;
        Bit_Err_mmse(i)=Bit_Err_mmse(i)+bit_err_mmse;
        Itr(i)=Itr(i)+1;
        disp("-------------");
        disp(Bit_Err_mmse./Itr./(n*k*(nOFDM-2)));
        disp(Bit_Err_zf./Itr./(n*k*(nOFDM-2)));
    end
    count=0;
    [BER_theor(i), SER_theor(i)]=berfading(EbNoArray(i),'qam',M,L);
end

BER_zf=Bit_Err_zf./Itr./(n*k*(nOFDM-2));
BER_mmse=Bit_Err_mmse./Itr./(n*k*(nOFDM-2));
SER=Sim_Err./Itr./(n*(nOFDM-2));
EVM_zf=evm_zf_avg./Itr;
EVM_mmse=evm_mmse_avg./Itr;


% Функция для вычисления EVM
function evm = calculate_evm(tx_qam, rx_fft_eq, n)
    % Вычисление EVM для ZF Equalizer
    evm = (1 / (n)) * sum((imag(tx_qam) - imag(rx_fft_eq)).^2 + (real(tx_qam) - real(rx_fft_eq)).^2, 'all');
end

function [qam_seq_bin, qam_list] = get_sequence(M)

    qam_seq=[sqrt(M)-1:-2:1 -1:-2:-(sqrt(M)-1)];
    qam_list=table2array(combinations(flip(qam_seq),qam_seq));

    qam_complex=zeros(M,1);
    for ii=1:length(qam_list)
        qam_complex(ii)=qam_list(ii,1)+i*qam_list(ii,2);
    end

    qam_seq_bin=qam_complex;
end

function [tx,qam_seq_gray] = qam_modulation(M,data,qam_seq,Es)
    symbols = bit2int(data,length(data));
    [~, seq_gray]=comm.internal.qam.getGrayMapping(M, symbols);
    qam_seq_gray=qam_seq(seq_gray+1);
    tx=qam_seq_gray(symbols+1)/sqrt(Es);
end

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