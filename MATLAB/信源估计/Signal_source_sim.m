B = 2*1e9; %带宽为2GHz
Fs = 6*1e9; %采样率6G

f0 = 2.5*1e5; %起始频率2.5*1e9;

T = 300;

[chirp_signal, t, N] = LFM_signal_gen(B, T, f0, 'R');
chirp_signal = sqrt(P)*chirp_signal;
figure(1);
subplot(2, 1, 1);
plot(t*1e6, real(chirp_signal));
xlabel('Time in u sec')
title('Real part of chirp signal');
grid on;
axis tight;
subplot(2,1,2);
freq=linspace(-fs/2, fs/2, N);
plot(freq*1e-6, fftshift(10*log10(abs(fft(chirp_signal)))));

