function singletone = singletoneGEN(fc,fs,power,N)
% 生成单频带干扰
%     参数设置：
%     fc：干扰中心功率
%     fs：采样率
%     power：功率
%     N：数据长度

%% ===================BPF 滤波器设计参数 ============================%

t=1/fs*(1:N);
ampJam  = sqrt(2*power);
singletone = ampJam*sin(2*pi*fc*t);
