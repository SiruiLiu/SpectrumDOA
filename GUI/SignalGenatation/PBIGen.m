function PBI = PBIGen(fc,Bw,fs,N)
% 生成部分频带干扰
%     参数设置：
%     fc：干扰中心功率
%     Bandwidth：干扰带宽
%     fs：采样率
%     power：功率
%     N：数据长度

%% ===================BPF 滤波器设计参数 ============================%

Wp = [fc-0.5*Bw  fc+0.5*Bw];       % 通带截至频率(Hz)
Ws = [fc-0.6*Bw  fc+0.6*Bw];       % 阻带截至频率(Hz)  过渡带指定为带宽的10%
Rp = 3;                            % 通带纹波(dB)
Rs = 80;                           % 阻带衰减(dB)

Wp = Wp/(fs/2);             % 频率归一
Ws = Ws/(fs/2);             % 频率归一

[ b, a ] = myremez( Wp, Ws, Rp, Rs );
jamFilter.Hd = dfilt.df1( b, a );                 % BPF滤波器

initialPhase = 0;                               % 载波初始相位
x = randn(N,1);    % 零均值高斯白噪声
% 干扰载波调制
s = repmat( x, 1, length( initialPhase ) );
theta = 2*pi*fc*[1:length(s)]'/fs;
theta = repmat( theta, 1, length( initialPhase ) ) + repmat( initialPhase.', length(s), 1 );
r = s.*cos( theta );

% 干扰成形滤波
jammerArray_temp1 = filter( jamFilter.Hd, r );

% 干扰功率归一
jammerArray_temp2 = jammerArray_temp1./repmat(std(jammerArray_temp1), length(s) ,1 );

% 设置干扰功率
PBI =  jammerArray_temp2';


return;
