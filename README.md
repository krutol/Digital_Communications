# Digital_Communications
In this project, a communication system consists of Bit source, transmitter, channel, receiver.

In the beginning, an input signal is generated as a sequence of bits by the Bit source. This sequence of bits is then transferred to transmitter from bit source. At the transmitter side, QPSK modulator is used to convert the sequence of bits into symbols. Then, the root raised cosine filter is used for pulse shaping with roll-off factor 0.35. After this, the output signal is sent through the channel.

In this project, we use two channels to investigation: Additive White Gaussian Noise (AWGN) and Rayleigh fading channel. The AWGN channel we use the MATLAB function 'awgn' to add to the channel. Rayleigh fading channel is modulated by multiplying the complex Gaussian random variable. At receiver side, Single and multiple antennas are employed for receiving the transmitted signal. Maximum ratio combining (MRC) which is a method of diversity combining, is used to combine the transmitted signals from single and multiple antennas.
