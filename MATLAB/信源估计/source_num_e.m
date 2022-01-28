%This program is used to estimate the source number
clc; clear all; close all;
%% define the parameter
lmd = 1; %set the wave length for 1m;
dis = 0.5 * lmd; % the array distance, this is the maximum for no aliasing
element_num = 8; % array element number
snapshot_num = 200; % snapshot number
snr = 0; % the SNR for per snapshot per array
doa = [5 30 90]; % direction of the signals
signal_num = length(doa); % signal number
fd = randi([1000 2000], 1, signal_num); % set the signal frequency
%fd = linspace(0, 2000, signal_num); % set the signal frequency
fs = 5000; % set the sample frequency
signal_power = 1; 
%signal power, all the signals power is equal, however,
%the error for the power can effort the resolution
step = 0.04; % step for search

%% array steering matrix
Asm = exp(1i * 2 * pi * (0 : element_num - 1)' * dis * cosd(doa) ./ lmd);
for snr = -20 : 5 : 30
    sp = zeros(1, 4);
    for i = 1 : 1000
        % generator the source signal
        rand_phase = 2 * pi * rand(1, signal_num); % signal random phase
        omega = 2 * pi * fd' ./ fs * (0 : snapshot_num - 1);
        % St = signal_power .^ 0.5 * exp(1i * omega) .* repmat(rand_phase', 1, snapshot_num);
        St = signal_power .^ 0.5 * exp(1i * omega);
        %generator the noise signal
        noise_power = signal_power ./ (10 .^ (snr / 10));
        noise = wgn(element_num, snapshot_num, noise_power, 'linear', 'complex');
        % the receive signal
        Xsr = Asm * St + noise;
        % covariance matrix
        Rx = Xsr * Xsr' ./ snapshot_num;

        % Information theory
        [V, D] = eig(Rx); %eigen decomposition
        Q = sort(diag(D), 'descend');% sort the eigenvalues in descending order
        % logarithm likelihood function
        L = zeros(1, element_num - 1);
        Paic = zeros(1, element_num - 1);
        Pmdl = zeros(1, element_num - 1);
        Phq = zeros(1, element_num - 1);
        for k = 1 : element_num - 1
            temp = element_num - k;
            Qk = Q(k + 1 : element_num);
            L(k) = snapshot_num * temp * log(1 / temp * sum(Qk) / power(prod(Qk), 1 / temp));
            Paic(k) = 2*k * (2 * element_num - k);
%             Pmdl(k) = k * (2 * element_num - k) * log (snapshot_num / 2);
            Pmdl(k) = k * (2 * element_num - k) * log (snapshot_num)*0.5;
            Phq(k) =  k * (2 * element_num - k) * log(log (snapshot_num / 2));
        end
        Kaic = L + Paic;
        Kmdl = L + Pmdl;
        Khq = L + Phq;
        kaic = find(Kaic == min(Kaic));
        if kaic == signal_num
            sp(1) = sp(1) + 1;
        end
        kmdl = find(Kmdl == min(Kmdl));
        if kmdl == signal_num
            sp(2) = sp(2) + 1;
        end
        khq = find(Khq == min(Khq));
        if khq == signal_num
            sp(3) = sp(3) + 1;
        end
        
        knew=SourceEst(Xsr);
        if knew == signal_num
            sp(4) = sp(4) + 1;
        end
    end
    n = snr / 5 + 5;
    s_p(n, :) = sp/10;
end
plot(-20 : 5 : 30, s_p(:, 1), '-*r')
hold on
plot(-20 : 5 : 30, s_p(:, 2), '-xb')
hold on
plot(-20 : 5 : 30, s_p(:, 3), '-dg')
plot(-20 : 5 : 30, s_p(:, 4), '-om')
hold off
xlabel('信噪比/dB')
ylabel('正确率/%');
title(['采样点数', num2str(snapshot_num)]);
legend('ALC', 'MDL', 'HQ', 'New');
grid on;