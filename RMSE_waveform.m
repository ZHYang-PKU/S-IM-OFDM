clc;
clear;
global c0 fc lambda M N delta_f Ts CPsize

%% 系统参数
c0 = 3e+8;  % 光速
fc = 5e+9; % 中心频点
lambda = c0 / fc; % 波长
M = 256; % 子载波总个数
K = 64;
N = 15; % 每子帧的符号个数
delta_f = 120e+3; % 子载波间隔
T = 1 / delta_f; % 符号周期
Tcp = T / 4;
Ts = T + Tcp; 
CPsize = M / 4;
bitsPerSymbol = 2; 
qam = 2^(bitsPerSymbol);
N_est=500; 

bitsPerSymbol_proposed=log2(K/16);

RMSE_r_OFDM=zeros(1,7);
RMSE_v_OFDM=zeros(1,7);
RMSE_r_proposed=zeros(1,7);
RMSE_v_proposed=zeros(1,7);


EbN0=-10:5:20;
tic
%% 
for i=1:1:7
    fprintf('EbN0(dB) is %g  \n',EbN0(i))
    error_r_OFDM=0;
    error_v_OFDM=0;
    error_r_SIM=0;
    error_v_SIM=0;
    for jj=1:1:N_est
        data = randi([0 qam - 1], M, N);
        TxData = qammod(data, qam, 'gray');

        % OFDM
        TxSignal_OFDM = ifft(TxData, M); % IFFT
        TxSignal_cp = [TxSignal_OFDM(M - CPsize + 1: M, :); TxSignal_OFDM]; % add CP
        TxSignal_cp = reshape(TxSignal_cp, [], 1); % time-domain transmit signal

        % SIMOFDM
        TxSignal_cp_SIM=SIM_Mod(data, qam, CPsize, M, K);

        %% target response
        SNR = EbN0(i);
        r = 60; 
        v = 40;  
        RxSignal_OFDM = sensingSignalGen(TxSignal_cp, r,v,SNR);
        RxSignal_SIM=sensingSignalGen_SIM(TxSignal_cp_SIM, r,v,SNR);
        k = length(r);

        %% 
        Rx = RxSignal_OFDM(1:size(TxSignal_cp,1),:); 
        Rx = reshape(Rx,[],N);
        Rx = Rx(CPsize + 1 : M + CPsize,:);  
        Rx_dem = fft(Rx,M);
        CIM_OFDM = Rx_dem .* conj(TxData);  
        [range_OFDM, velocity_OFDM]=OFDMsensing(CIM_OFDM,k, i);
        error_r_OFDM=error_r_OFDM +(range_OFDM-r)^2;
        error_v_OFDM=error_v_OFDM + (velocity_OFDM-v)^2;

        %% 
        Rx = RxSignal_SIM(1:size(TxSignal_cp,1),:); 
        Rx = reshape(Rx,[],N);
        Rx = Rx(CPsize + 1 : M + CPsize,:); 
        Rx_dem = fft(Rx,M);
        
        CIM = Rx_dem .*conj(TxData); 
        [range,velocity] = SIMsensing(CIM,k,i);
        error_r_SIM=error_r_SIM + (range-r)^2;
        error_v_SIM=error_v_SIM + (velocity-v)^2;
    end
    RMSE_r_OFDM(i)=sqrt(error_r_OFDM/(N_est));
    RMSE_v_OFDM(i)=sqrt(error_v_OFDM/(N_est));
    RMSE_r_proposed(i)=sqrt(error_r_SIM/(N_est));
    RMSE_v_proposed(i)=sqrt(error_v_SIM/(N_est));
    
end
toc
figure(1);
semilogy(EbN0,RMSE_r_proposed, 'rp-','LineWidth',1.8,'MarkerSize',8 );hold on;
semilogy(EbN0,RMSE_r_OFDM, 'bo-','LineWidth',1.8,'MarkerSize',8 ); hold on;
grid on; box on;
xlabel('SNR(dB)');
ylabel('RMSE-range (m)');
legend('S-IM-OFDM','OFDM');
ylim([0.001,10]);

figure(2);
semilogy(EbN0,RMSE_v_proposed, 'rp-','LineWidth',1.8,'MarkerSize',8 );hold on;
semilogy(EbN0,RMSE_v_OFDM, 'bo-','LineWidth',1.8,'MarkerSize',8 ); hold on;
grid on; box on;
xlabel('SNR(dB)');
ylabel('RMSE-velocity (m)');
legend('S-IM-OFDM','OFDM');
ylim([0.001,10]);
