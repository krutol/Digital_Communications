clear 
clc
close all



M = 4;
Ph_err = pi/6;
k = log2(M);


EbNoVec = -10:2:30;
berAWGN = zeros(size(EbNoVec));
berRay = zeros(size(EbNoVec));
numSum = 1e5;

%% AWGN Phase est

for i = 1:length(EbNoVec)
    
    snrdb = EbNoVec(i) + 10*log10(k);

        dataIn = randi([0,M-1],numSum,k);

        txSig0 = pskmod(dataIn, M,Ph_err);
        txSig0_noPhase = pskmod(dataIn, M);
        rxSig0 = awgn(txSig0,snrdb,'measured');
        rxSig0_noPhase = awgn(txSig0_noPhase,snrdb,'measured');
        
        % phase est
        x4 = pskmod([0,M-1],M);
        x4_conj = mean((conj(x4)).^4);
        Rx_4 = mean((rxSig0).^4);
        Ph_err_est = 0.25*mean(angle(x4_conj*Rx_4));
        
        
        
        estArray(i) = Ph_err_est;
        
        rxSym0 = pskdemod(rxSig0_noPhase, M);
        
        rxSym0_beforeEst = pskdemod(rxSig0, M);
        
        rxSym0_est = pskdemod(rxSig0, M,Ph_err_est);
        
        [numErr0,berAWGN(i)] = biterr(dataIn, rxSym0_beforeEst);
        
        [numEr,berAWGN_noPhase(i)] = biterr(dataIn,rxSym0);
        
        [numE, berAWGN_est(i)] = biterr(dataIn,rxSym0_est);
end


%% Rayleigh Phase est

for i = 1: length(EbNoVec)
    
    snrdb = 10.^(0.1*EbNoVec(i));

        dataIn0 = randi([0,M-1], numSum, k);
 
        txSig1 = pskmod(dataIn0, M,Ph_err);
        noiseRay = 1/sqrt(2*snrdb*k)*(randn(length(txSig1),1) + 1i*randn(length(txSig1),1));% Complex Gaussian niose

        Rayleigh = 1/sqrt(2)*(randn(length(txSig1),1) + 1i*randn(length(txSig1),1));%Rayleigh

        rxSig1 = (Rayleigh.*txSig1+ noiseRay)./Rayleigh;%Rayleigh fading : output of channel
        
        % phase est
        x4 = pskmod([0,M-1],M);
        x4_conj = mean((conj(x4)).^4);
        Rx_4 = mean((rxSig1).^4);
        Ph_err_est = abs(0.25*mean(angle(x4_conj*Rx_4)));
        
        estArrayRay(i) = abs(Ph_err_est);
        
        
        rxSym1 = pskdemod(rxSig1 , M,Ph_err);
        
        rxSym1_beforeEst = pskdemod(rxSig1, M);
        
        rxSym1_est = pskdemod(rxSig1, M,Ph_err_est);

        
        [nError1,berRay(i)] = biterr(dataIn0, rxSym1);
        [nError2, berRay_beforeEst(i)] = biterr(dataIn0, rxSym1_beforeEst);
        [nError3, berRay_Est(i)] = biterr(dataIn0, rxSym1_est);  

    
end









%% 
figure(1)
semilogy(EbNoVec,berAWGN_noPhase, EbNoVec,berAWGN_est,':X',EbNoVec,berAWGN)
axis([-10,20,1e-5,0.5]);
grid on;
legend('AWGNnoPhase','AWGNPhaseEst','AWGNbeforeEst');
xlabel('Eb/No(dB)');
ylabel('BER');
title('BER for AWGN');



figure(2)


semilogy(EbNoVec,berRay, EbNoVec,berRay_beforeEst,EbNoVec,berRay_Est,':rX')
axis([-10,20,1e-5,0.5]);
grid on;
legend('Ray','RaybeforeEst','RayPhaseEst');
xlabel('Eb/No(dB)');
ylabel('BER');
title('BER for Ray');

figure(3)
plot(EbNoVec,estArray)
hold on;
plot(EbNoVec,estArrayRay)
legend('AWGN-estError','Ray-estError');
xlabel('Eb/No(dB)');
ylabel('Phase estimate');
title('Estimated Phase based on maximum-likelihood in both channels');


% figure(4)
% plot(abs(estArrayRay-Ph_err))



