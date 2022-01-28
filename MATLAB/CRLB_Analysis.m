clear;
close all;
clc;

N = 8;
r=0.15;
C = 299792458;
Freq=1575e6;

Lambda = C/Freq;
TABLE = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 16, 24, 32, 64, 128, 256, 512, 10^3, 1024, 2048, 4096, 8192, 10^4];
SNR0 = zeros(1, length(TABLE));
SNR1 = zeros(1, length(TABLE));
SNR3 = zeros(1, length(TABLE));
SNR2 = zeros(1, length(TABLE));

for i = 1:length(SNR0)
    SNR0(i) = calCRLB(10, TABLE(i), N, r, Lambda);
    SNR1(i) = calCRLB(1, TABLE(i), N, r, Lambda);
    SNR2(i) = calCRLB(0.1, TABLE(i), N, r, Lambda);
    SNR3(i) = calCRLB(0.01, TABLE(i), N, r, Lambda);
end

SNR0_1 = SNR0*180/pi;
SNR1_1 = SNR1*180/pi;
SNR2_1 = SNR2*180/pi;
SNR3_1 = SNR3*180/pi;

figure; 

x = logspace(0, 4, 1000);
y1 = calCRLB(10, x, N, r, Lambda);
y2 = calCRLB(1, x, N, r, Lambda);
y3 = calCRLB(0.1, x, N, r, Lambda);
y4 = calCRLB(0.01, x, N, r, Lambda);
loglog(x, y1, x, y2, x, y3, x, y4);
yticks([0, 1, 5, 10, 30, 60]);
grid on;
legend('SNR=10dB', 'SNR=0dB', 'SNR=-10dB', 'SNR=-20dB');
ylabel('误差/度');
xlabel('快拍数');