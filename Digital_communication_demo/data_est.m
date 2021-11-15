clear 
clc
close all

M = 4;
EbNoVec = 15;
numSum = 100:50:10000;
ph_err = pi/6;
k = log2(M);


%% 2 antenna%%%%%%%


for j = 1 : length(numSum)
    snrdb = EbNoVec+10*log(k);
    dataIn = randi([0,M-1], numSum(j) , k);
    txRig = pskmod(dataIn, M,ph_err,'gray');
    
 
    rxSig_a1 = awgn(txRig,snrdb,'measured');
    rxSig_a2 =  awgn(txRig,snrdb,'measured');
    rxSig_tot2 = 1/2*(rxSig_a1 + rxSig_a2);
    
    
    % phase est
        x4 = pskmod([0,M-1],M);
        x4_conj = mean((conj(x4)).^4);
        Rx_4 = mean((rxSig_a1).^4);
        Ph_err_est_a1 = 0.25*mean(angle(x4_conj*Rx_4));
        
     % phase est
        x5 = pskmod([0,M-1],M);
        x5_conj = mean((conj(x5)).^4);
        Rx_5 = mean((rxSig_tot2).^4);
        Ph_err_est_a1a2 = 0.25*mean(angle(x5_conj*Rx_5));
        
        
     Sub_01(j) = abs(Ph_err_est_a1 - ph_err);
     
     Sub_02(j) = abs(Ph_err_est_a1a2 - ph_err);
     

    
    
    
    
    
    
    
    
end

%% ray_2antenna

for i = 1 : length(numSum)
    snrdb = 10.^(0.1*EbNoVec);
    dataInRay2 = randi([0,M-1], numSum(i) , k);
    txRigRay2 = pskmod(dataInRay2, M,ph_err,'gray');
    
    %%%% #1
    noiseRay_a1 = 1/sqrt(2*snrdb*k)*(randn(length(txRigRay2),1) + 1i*randn(length(txRigRay2),1));% Complex Gaussian niose

    Rayleigh_a1 = 1/sqrt(2)*(randn(length(txRigRay2),1) + 1i*randn(length(txRigRay2),1));%Rayleigh
    
    Rayleigh_nominator_a1 = Rayleigh_a1.*txRigRay2+ noiseRay_a1;

    rxSigRay_one = Rayleigh_nominator_a1./Rayleigh_a1;%Rayleigh fading : output of channel
    
    %%%% #2
    
    noiseRay_a2 = 1/sqrt(2*snrdb*k)*(randn(length(txRigRay2),1) + 1i*randn(length(txRigRay2),1));% Complex Gaussian niose

    Rayleigh_a2 = 1/2*1/sqrt(2)*(randn(length(txRigRay2),1) + 1i*randn(length(txRigRay2),1));%Rayleigh
    
    Rayleigh_nominator_a2 = Rayleigh_a2.*txRigRay2+ noiseRay_a2;

    rxSigRay_two = Rayleigh_nominator_a2./Rayleigh_a2;%Rayleigh fading : output of channel
    
    %%%% combine
    
    nominator = conj(Rayleigh_a1) .* Rayleigh_nominator_a1 + conj(Rayleigh_a2) .* Rayleigh_nominator_a2;
    
    denominator = conj(Rayleigh_a1) .* Rayleigh_a1 + conj(Rayleigh_a2) .* Rayleigh_a2;
    
    equalizer = nominator./denominator;
    
     % phase est
        x4 = pskmod([0,M-1],M);
        x4_conj = mean((conj(x4)).^4);
        Rx_4 = mean((rxSigRay_one).^4);
        Ph_err_RayEst_a1 = abs(0.25*mean(angle(x4_conj*Rx_4)));
        
        
     % phase est
        x5 = pskmod([0,M-1],M);
        x5_conj = mean((conj(x5)).^4);
        Rx_5 = mean((equalizer).^4);
        Ph_err_RayEst_a1a2 = 0.25*mean(angle(x5_conj*Rx_5));
    
     Sub_Ray01(i) = abs(Ph_err_RayEst_a1 - ph_err);
     
     Sub_Ray02(i) = abs(Ph_err_RayEst_a1a2 - ph_err);
    
    
end


%% AWGN 4 antenna

for j = 1 : length(numSum)
    snrdb = EbNoVec + 10*log(k);
    dataIn = randi([0,M-1], numSum(j) , k);
    txRigAWGN4 = pskmod(dataIn, M,ph_err,'gray');
    
  
    rxSigAWGN_a1 = awgn(txRigAWGN4,snrdb,'measured');
    rxSigAWGN_a2 = awgn(txRigAWGN4,snrdb,'measured');
    rxSigAWGN_a3 = awgn(txRigAWGN4,snrdb,'measured');
    rxSigAWGN_a4 = awgn(txRigAWGN4,snrdb,'measured');
    
    rxSigAWGN_tot4 = 1/4*(rxSigAWGN_a1 + rxSigAWGN_a2 + rxSigAWGN_a3 + rxSigAWGN_a4);
    
    
    % phase est
        x5 = pskmod([0,M-1],M);
        x5_conj = mean((conj(x5)).^4);
        Rx_5 = mean((rxSigAWGN_tot4).^4);
        Ph_err_est4 = 0.25*mean(angle(x5_conj*Rx_5));
    
        
        Sub_04(j) = abs(Ph_err_est4 - ph_err);
    
    
    
    
    
end





%% Ray 4 antenna


for i = 1 : length(numSum)
    snrdb = 10.^(0.1*EbNoVec);
    dataInRay4 = randi([0,M-1], numSum(i) , k);
    txRigRay4 = pskmod(dataInRay4, M,ph_err,'gray');
    
    %%%% #1
    noiseRay_a1 = 1/sqrt(2*snrdb*k)*(randn(length(txRigRay4),1) + 1i*randn(length(txRigRay4),1));% Complex Gaussian niose

    Rayleigh_a1 = 1/sqrt(2)*(randn(length(txRigRay4),1) + 1i*randn(length(txRigRay4),1));%Rayleigh
    
    Rayleigh_nominator_a1 = Rayleigh_a1.*txRigRay4+ noiseRay_a1;

    rxSigRay_one = Rayleigh_nominator_a1./Rayleigh_a1;%Rayleigh fading : output of channel
    
    %%%% #2
    
    noiseRay_a2 = 1/sqrt(2*snrdb*k)*(randn(length(txRigRay4),1) + 1i*randn(length(txRigRay4),1));% Complex Gaussian niose

    Rayleigh_a2 = 1/2*1/sqrt(2)*(randn(length(txRigRay4),1) + 1i*randn(length(txRigRay4),1));%Rayleigh
    
    Rayleigh_nominator_a2 = Rayleigh_a2.*txRigRay4+ noiseRay_a2;

    rxSigRay_two = Rayleigh_nominator_a2./Rayleigh_a2;%Rayleigh fading : output of channel
    
    
    %%%% #3
    
    noiseRay_a3 = 1/sqrt(2*snrdb*k)*(randn(length(txRigRay4),1) + 1i*randn(length(txRigRay4),1));% Complex Gaussian niose

    Rayleigh_a3 = 1/3*1/sqrt(2)*(randn(length(txRigRay4),1) + 1i*randn(length(txRigRay4),1));%Rayleigh
    
    Rayleigh_nominator_a3 = Rayleigh_a3.*txRigRay4+ noiseRay_a3;

    rxSigRay_three = Rayleigh_nominator_a3./Rayleigh_a3;%Rayleigh fading : output of channel
    
    
    %%%% #4
    
    noiseRay_a4 = 1/sqrt(2*snrdb*k)*(randn(length(txRigRay4),1) + 1i*randn(length(txRigRay4),1));% Complex Gaussian niose

    Rayleigh_a4 = 1/4*1/sqrt(2)*(randn(length(txRigRay4),1) + 1i*randn(length(txRigRay4),1));%Rayleigh
    
    Rayleigh_nominator_a4 = Rayleigh_a4.*txRigRay4+ noiseRay_a4;

    rxSigRay_four = Rayleigh_nominator_a4./Rayleigh_a4;%Rayleigh fading : output of channel
    
    %%%% combine
    
    nominator = conj(Rayleigh_a1) .* Rayleigh_nominator_a1 + conj(Rayleigh_a2) .* Rayleigh_nominator_a2 + conj(Rayleigh_a3) .* Rayleigh_nominator_a3 + conj(Rayleigh_a4) .* Rayleigh_nominator_a4;
    
    denominator = conj(Rayleigh_a1) .* Rayleigh_a1 + conj(Rayleigh_a2) .* Rayleigh_a2 + conj(Rayleigh_a3) .* Rayleigh_a3 + conj(Rayleigh_a4) .* Rayleigh_a4;
    
    equalizer = nominator./denominator;
    
    % phase est
        x5 = pskmod([0,M-1],M);
        x5_conj = mean((conj(x5)).^4);
        Rx_5 = mean((equalizer).^4);
        Ph_err_RayEst4 = 0.25*mean(angle(x5_conj*Rx_5));
        
        
   Sub_Ray04(i) = abs(Ph_err_RayEst4 - ph_err);
end

figure(1)
plot(numSum,Sub_01,numSum,Sub_02,numSum,Sub_04)
grid on;
legend('1 antenna','2 antenna','4 antenna');
xlabel('numSum');
ylabel('phase error');
axis([100,10000,0,0.005]);
title('AWGN Inphase error for diversity antenna');

figure(2)

plot(numSum,Sub_Ray01,numSum,Sub_Ray02,numSum,Sub_Ray04)
grid on;
legend('1 antenna','2 antenna','4 antenna');
xlabel('numSum');
ylabel('phase error');
axis([100,10000,0,0.05]);
title('Rayleigh Inphase error for diversity antenna');





