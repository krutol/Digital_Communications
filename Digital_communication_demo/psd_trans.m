clear 
close all


%Tong jiang, ziyang li, shuo zhang
M = 4;
k = log2(M);
Ph_err = pi/4;
symbolrate = 7e3;
rolloff = 0:0.35:1;
data = randi([0 M-1], 10000, 1);
mod = pskmod(data, M, Ph_err);
span = 6;%6 symbols
sps = 4;%4samples
for i= 1:length(rolloff)
    rcosfilter = rcosdesign(rolloff(i),span,sps,'sqrt');
    upsample = upfirdn(mod, rcosfilter, sps);%upsampling
    
    figure(1);
    pwelch(upsample,hamming(1024),[],[],symbolrate*sps,'centered');
    hold on;
    grid on;
    legend('Roll Off : 0','Roll Off : 0.35','Roll Off : 0.7');
    xlabel('Freq');
    ylabel('Power');
    title('PSD FOR QPSK');
end






%% CONSTELLATION%%%%%%%%%


rxSig = awgn(mod,20);
figure(2)
scatterplot(rxSig)
grid on
title('signal constellation in AWGN when SNR = 20 dB')%constellation in 20db


rxSig1 = awgn(mod,10);
figure(3)
scatterplot(rxSig1)
grid on
title('signal constellation in AWGN when SNR = 10 dB')%constellation in -10db


rxSig2 = awgn(mod,0);
figure(4)
scatterplot(rxSig2)
grid on
title('signal constellation in AWGN when SNR = 0 dB')%constellation in 0db


rxSig3 = awgn(mod,-10);
figure(5)
scatterplot(rxSig3)
grid on
title('signal constellation in AWGN when SNR = -10 dB')%constellation in 10db


%% 0ne antenna compare theo AWGN and Ray%%%%%%%%

EbNoVec = -10:2:30;
berAWGN = zeros(size(EbNoVec));
berRay = zeros(size(EbNoVec));
numSum = 1e5;


%% ber awgn%%%%
for i = 1:length(EbNoVec)
    
    snrdb = EbNoVec(i) + 10*log10(k);

        dataIn = randi([0,M-1],numSum,k);

        txSig0 = pskmod(dataIn, M, Ph_err,'gray');

        rxSig0 = awgn(txSig0,snrdb,'measured');
        
        rxSym0 = pskdemod(rxSig0, M, Ph_err, 'gray');

        [numErr0,berAWGN(i)] = biterr(dataIn, rxSym0);
        

end
    

     



%% ber ray%%%

for i = 1: length(EbNoVec)
    
    snrdb = 10.^(0.1*EbNoVec(i));

        dataIn0 = randi([0,M-1], numSum, k);
 
        txSig1 = pskmod(dataIn0, M,Ph_err,'gray');
        noiseRay = 1/sqrt(2*snrdb*k)*(randn(length(txSig1),1) + 1i*randn(length(txSig1),1));% Complex Gaussian niose

        Rayleigh = 1/sqrt(2)*(randn(length(txSig1),1) + 1i*randn(length(txSig1),1));%Rayleigh

        rxSig1 = (Rayleigh.*txSig1+ noiseRay)./Rayleigh;%Rayleigh fading : output of channel
        rxSym1 = pskdemod(rxSig1, M,Ph_err,'gray');

        
        [nError1,berRay(i)] = biterr(dataIn0, rxSym1);

    
end






%% theo%%%%
theoAWGN = berawgn(EbNoVec, 'psk', M,'nondiff');
theoRay = berfading(EbNoVec, 'psk',M , 1);%1 is diversity order
figure(7)
semilogy(EbNoVec,theoAWGN,EbNoVec, berAWGN,'*',EbNoVec,theoRay,EbNoVec, berRay,'*')
axis([-10,20,1e-5,0.5]);
grid on;
legend('TheoryAWGN','newAWGN','TheoryRay','newRay');
xlabel('Eb/No(dB)');
ylabel('BER');
title('BER for theory AWGN and Ray');



