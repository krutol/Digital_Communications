clear 
clc
close all

M = 4;
EbNoVec = -10:2:30;
numSum = 1e5;
ph_err = pi/6;
k = log2(M);


%% 2 antenna%%%%%%%

BerAWGN_2Ant = [];
BerAWGN_1Ant = [];
BerRay_2Ant = [];
BerRay_1Ant = [];


for j = 1 : length(EbNoVec)
    snrdb = EbNoVec(j)+10*log(k);
    dataIn = randi([0,M-1], numSum , k);
    txRig = pskmod(dataIn, M,ph_err,'gray');
    
 
    rxSig_a1 = awgn(txRig,snrdb,'measured');
    rxSig_a2 =  awgn(txRig,snrdb,'measured');
    rxSig_tot2 = 1/2*(rxSig_a1 + rxSig_a2);
    
    rxSym_a1 = pskdemod(rxSig_a1,M,ph_err,'gray');
    rxSym_tot = pskdemod(rxSig_tot2,M,ph_err,'gray');
    
    [n1,BerAWGN_1Ant(j)] = biterr(rxSym_a1,dataIn);
    [n2,BerAWGN_2Ant(j)] = biterr(rxSym_tot, dataIn);
    
    
    
    
    
    
    
    
end

%% ray_2antenna

for i = 1 : length(EbNoVec)
    snrdb = 10.^(0.1*EbNoVec(i));
    dataInRay2 = randi([0,M-1], numSum , k);
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
    
    demodDataRay2 = pskdemod(equalizer,M,ph_err,'gray');
    
    demodDataRay1 = pskdemod(rxSigRay_one,M,ph_err,'gray');
    
     [nError1,BerRay_2Ant(i)] = biterr(dataInRay2, demodDataRay2);
     
     [nError2,BerRay_1Ant(i)] = biterr(dataInRay2, demodDataRay1);
    
    
    
    
end


%% AWGN 4 antenna
BerAWGN_4Ant = [];
for j = 1 : length(EbNoVec)
    snrdb = EbNoVec(j)+10*log(k);
    dataIn = randi([0,M-1], numSum , k);
    txRigAWGN4 = pskmod(dataIn, M,ph_err,'gray');
    
  
    rxSigAWGN_a1 = awgn(txRigAWGN4,snrdb,'measured');
    rxSigAWGN_a2 = awgn(txRigAWGN4,snrdb,'measured');
    rxSigAWGN_a3 = awgn(txRigAWGN4,snrdb,'measured');
    rxSigAWGN_a4 = awgn(txRigAWGN4,snrdb,'measured');
    
    rxSigAWGN_tot4 = 1/4*(rxSigAWGN_a1 + rxSigAWGN_a2 + rxSigAWGN_a3 + rxSigAWGN_a4);
    
    
    rxSymAWGN_tot4 = pskdemod(rxSigAWGN_tot4,M,ph_err,'gray');
    
    [n2,BerAWGN_4Ant(j)] = biterr(rxSymAWGN_tot4, dataIn);
    
    
    
    
    
    
end





%% Ray 4 antenna

BerRay_4Ant = [];
for i = 1 : length(EbNoVec)
    snrdb = 10.^(0.1*EbNoVec(i));
    dataInRay4 = randi([0,M-1], numSum , k);
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
    
    demodDataRay4 = pskdemod(equalizer,M,ph_err,'gray');
    
 
    
     [nError1,BerRay_4Ant(i)] = biterr(dataInRay4, demodDataRay4);
   
end


theoRay1 = berfading(EbNoVec, 'psk',M , 1);%1 is diversity order
theoRay2 = berfading(EbNoVec, 'psk',M , 2);%2 is diversity order
theoRay4 = berfading(EbNoVec, 'psk',M , 4);%4 is diversity order







figure(1)
semilogy(EbNoVec,BerRay_2Ant,'r',EbNoVec,theoRay2,':rX',EbNoVec,BerRay_1Ant,'b',EbNoVec,theoRay1,':bX',EbNoVec,BerRay_4Ant,'k',EbNoVec,theoRay4,':kX')
grid on;
legend('2 antenna','theo2','1 antenna','theo1','4 antenna','theo4');
xlabel('Eb/No');
ylabel('BER');
axis([-10,20,1e-5,0.5]);
title('Rayleigh BER for diversity antenna');

figure(2)

semilogy(EbNoVec,BerAWGN_1Ant,EbNoVec,BerAWGN_2Ant,EbNoVec,BerAWGN_4Ant);
grid on;
legend('1 antenna','2 antenna','4 antenna');
xlabel('Eb/No');
ylabel('BER');
axis([-10,20,1e-5,0.5]);
title('AWGN BER for diversity antenna');


