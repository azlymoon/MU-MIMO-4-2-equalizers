% Сохраните этот код в файл с расширением .m, например, main.m

clear all
fe=12000;       %the number of collected statistics at each noise var value

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
rb = 100;                % Number of RB (resource blocks)
n = rb*12;               % Number of QAM symbols in OFDM symbol
padding_len=(nFFT-n)/2;  % length of padding zeros elements before ifft (half of them)
CP_len=0.125*nFFT;       % length of cyclic prefix (guard band)

%Fading channel parameters
fs = 30.72e6;            % Sample rate in Hz from LTE

%Number of antennas
L = 1;

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

%chs.PRBSet = (0:4).';
chs.PRBSet = (0:49).';

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

% Параметры временного фильтра (DFE) - как указано в вашем коде
L_kih = 40;

% Параметры блока C|P (Comparator/Processor).  Простое пороговое устройство.
threshold = 0;

%Converting EsNO from EbNo
EsNoArray=EbNoArray+10*log10(k);%+10*log10(n/nFFT);

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
        %FFT step, without DM-RS
        norm_1 = 1/sqrt(n);
        tx_qam_shift = fftshift(tx_qam,1);
        tx_qam_fft = norm_1.*fft(tx_qam_shift,n,1);

        %Add dm-rs signal on positions sym_1 and sym_2
        tx_qam_dmrs=zeros(n,nOFDM);

        tx_qam_dmrs(:,[1:sym_1-1,sym_1+1:sym_2-1,sym_2+1:end])=tx_qam_fft;
        tx_qam_dmrs(:,sym_1)=dmrs_1;
        tx_qam_dmrs(:,sym_2)=dmrs_2;

        dmrs_fft_1 = norm_1.*fft(dmrs_1);
        dmrs_fft_2 = norm_1.*fft(dmrs_2);

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

        %%%%%AWGN channel%%%%%%
        rx = zeros(size(rx_fading_lte));
        for l = 1:L
            noise = (randn(size(rx_fading_lte(:, l))) + 1i * randn(size(rx_fading_lte(:, l)))) / sqrt(2);
            noise((nFFT + CP_len) * (sym_1 - 1) + 1 + delay:(nFFT + CP_len) * sym_1 + delay) = complex(zeros(1, nFFT + CP_len));
            noise((nFFT + CP_len) * (sym_2 - 1) + 1 + delay:(nFFT + CP_len) * sym_2 + delay) = complex(zeros(1, nFFT + CP_len));
            rx(:, l) = rx_fading_lte(:, l) + sqrt(noise_var) * noise;
        end

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
        W_zf = zeros([n, L, 2]);
        for jj = 1:2
            for ii = 1:n
                W_zf(ii, :, jj) = pinv(h_channel(:, ii, jj));
            end
        end

        rx_fft_eq_zf=zeros(n,nOFDM);
        for j=1:7
            for ii=1:n
                rx_fft_eq_zf(ii,j)=W_zf(ii,:,1)*squeeze(rx_fft(ii,j,:));
                rx_fft_eq_zf(ii,j+7)=W_zf(ii,:,2)*squeeze(rx_fft(ii,j+7,:));
            end
        end

        %MMSE equalizer
        W_mmse = zeros([n, L, 2]);
        for jj = 1:2
            for ii = 1:n
                W_mmse(ii, :, jj) = ((h_channel(:, ii, jj)' * h_channel(:, ii, jj) + noise_var))\h_channel(:, ii, jj)';
            end
        end

        rx_fft_eq_mmse=zeros(n,nOFDM);

        %for j=1:7
        %    for ii=1:n
        %        rx_fft_eq_mmse(ii,j)=W_mmse(ii,:,1)*squeeze(rx_fft(ii,j,:));
        %        rx_fft_eq_mmse(ii,j+7)=W_mmse(ii,:,2)*squeeze(rx_fft(ii,j+7,:));
        %    end
        %end

        % Pre-allocate g_fb_matrix and G_FF_matrix with the correct size
        g_fb_matrix = zeros(L_kih+1, 2, L); % L_kih rows, 2 columns (for each DM-RS), L layers
        G_FF_matrix = zeros(n, 2, L);

        rx_fft_eq_zf_shift=ifft(rx_fft_eq_zf(:,[1:sym_1-1,sym_1+1:sym_2-1,sym_2+1:end]),n,1)./norm_1;
        rx_fft_eq_zf = ifftshift(rx_fft_eq_zf_shift,1);
        for l = 1:L  %Итерируемся по антеннам для расчета DFE коэффициентов
            for jj = 1:2 %итерируемся по номерам DM-RS
                [g_fb, G_FF]=Coeffients_DFE_calculation(squeeze(W_mmse(:,l,jj)),h_channel(l, :, jj), noise_var, L_kih, n);
                g_fb = [0+1i*0; g_fb];
                g_fb_matrix(:,jj,l) = g_fb; % сохраняем коэф. для каждой антенны
                G_FF_matrix(:,jj,l) = G_FF;
            end
        end
        for j=1:7
            for ii=1:n
                 rx_fft_eq_mmse(ii,j)=squeeze(G_FF_matrix(ii,1,:)).'*squeeze(rx_fft(ii,j,:));
                 rx_fft_eq_mmse(ii,j+7)=squeeze(G_FF_matrix(ii,2,:)).'*squeeze(rx_fft(ii,j+7,:));
            end
        end

        rx_fft_tde_shift=ifft(rx_fft_eq_mmse(:,[1:sym_1-1,sym_1+1:sym_2-1,sym_2+1:end]),n,1)./norm_1;
        rx_fft_tde = ifftshift(rx_fft_tde_shift,1);
        [~,rx_fft_dfe] = FIR_filter(rx_fft_tde,norm_1,n,nOFDM,g_fb_matrix,L_kih,M,sym_1,sym_2);
        [rx_fft_dfe_demod,rx_fft_dfe] = FIR_filter(rx_fft_dfe,norm_1,n,nOFDM,g_fb_matrix,L_kih,M,sym_1,sym_2);

        evm_zf=zeros(1,nOFDM-2);
        evm_mmse=zeros(1,nOFDM-2);
        for jj=1:nOFDM-2
            evm_zf(jj)=calculate_evm(tx_qam(:,jj),rx_fft_tde(:,jj),n);
            evm_mmse(jj)=calculate_evm(tx_qam(:,jj),rx_fft_tde(:,jj),n,h_channel(1,:,1));
        end

        evm_zf_avg(i)=evm_zf_avg(i)+mean(10*log10(evm_zf));
        evm_mmse_avg(i)=evm_mmse_avg(i)+mean(10*log10(evm_mmse));

        dw_zf=zeros(n*k,nOFDM-2);
        dw_const=zeros(n,nOFDM-2);
        dw_mmse=zeros(n*k,nOFDM-2);
        for jj=1:nOFDM-2
            for ii=0:n-1
                dw_mmse(ii*k+1:ii*k+k,jj)=qamdemod(rx_fft_dfe_demod(ii+1,jj),M,'gray',OutputType='bit',UnitAveragePower=true);
                dw_zf(ii*k+1:ii*k+k,jj)=qamdemod(rx_fft_eq_zf(ii+1,jj),M,'gray',OutputType='bit',UnitAveragePower=true);
                %qam_point = qamdemod(rx_fft_eq_zf(ii+1,jj),M,'gray',OutputType='approxllr',UnitAveragePower=true);
                %dw_const(ii+1,jj) = qam_point(1)+1i*qam_point(2);
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
    end
    disp("-------------");
    disp(Bit_Err_mmse./Itr./(n*k*(nOFDM-2)));
    disp(Bit_Err_zf./Itr./(n*k*(nOFDM-2)));
    count=0;
    [BER_theor(i), SER_theor(i)]=berfading(EbNoArray(i),'qam',M,L);
end

BER_zf=Bit_Err_zf./Itr./(n*k*(nOFDM-2));
BER_mmse=Bit_Err_mmse./Itr./(n*k*(nOFDM-2));
SER=Sim_Err./Itr./(n*(nOFDM-2));
EVM_zf=evm_zf_avg./Itr;
EVM_mmse=evm_mmse_avg./Itr;

% Функция для вычисления EVM
function evm = calculate_evm(tx_qam, rx_fft_eq, n, h_channel)
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

function [g_fb, G_FF]=Coeffients_DFE_calculation(W,h_channel, noise_var, L_kih, n)
    b_mmse = zeros(1,L_kih);
    M=length(h_channel);
    for j=1:L_kih
        b_curr = 0;
        for k=1:M
            b_curr = b_curr + exp((2*pi*1i*j*k)/M)/(abs(h_channel(k))^2+noise_var);
        end
        b_mmse(j) = -b_curr;
    end

    A_mmse = zeros(L_kih,L_kih);
    for j=1:L_kih
        for l=1:L_kih
            A_curr = 0;
            for k=1:M
                A_curr = A_curr + exp((2*pi*1i*(l-j)*k)/M)/(abs(h_channel(k))^2+noise_var);
            end
            A_mmse(j,l) = - A_curr;
        end
    end

    g_fb = linsolve(A_mmse,b_mmse.');
    g_fb_padding = [zeros(n/2-L_kih/2,1); g_fb; zeros(n/2-L_kih/2,1)];
    G_FB = (1/sqrt(n)).*fft(g_fb_padding);
    G_FF = W.*(1+G_FB);
end

function [rx_fft_dfe_demod,rx_fft_dfe] = FIR_filter(rx_fft_tde,norm_1,n,nOFDM,g_fb_matrix,L_kih,M,sym_1,sym_2)
    
    rx_fft_dfe = zeros(size(rx_fft_tde));
    rx_fft_dfe_demod = zeros(size(rx_fft_tde));
    
    for jj=1:nOFDM-2
        for ii=1:n
            if jj < 7
                g_fb = g_fb_matrix(:,1,1);
            else
                g_fb = g_fb_matrix(:,2,1);
            end
    
            symbol_window = zeros(1, L_kih+1);
            for k_tde = 1:L_kih+1
                if (ii-k_tde)<0
                    break
                end
                symbol_window(k_tde) = rx_fft_dfe_demod(ii-k_tde+1,jj);
            end
            
            % Свертка с g_fb (применение DFE)
            feedback_signal = sum(g_fb.' .* symbol_window);
    
            %Вычитание обратной связи
            rx_fft_dfe(ii,jj) = rx_fft_tde(ii,jj) - feedback_signal;
            qam_point = qamdemod(rx_fft_dfe(ii,jj),M,'gray',OutputType='bit',UnitAveragePower=true);
            rx_fft_dfe_demod(ii,jj)=qammod(qam_point,M,'gray',InputType='bit',UnitAveragePower=true);
        end
    end
end