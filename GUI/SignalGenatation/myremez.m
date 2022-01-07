function [ b, a ] = myremez( Wp, Ws, Rp, Rs )
% 采用parks-mcclellan方法计算FIR滤波器的阶数，并用Remez算法进行最优FIR滤波器设计
% 输入参数说明：
% Wp                    通带截至频率，如[ 0.3, 0.7 ]，定义方法与matlab相同
% Ws                    阻带截至频率，如[ 0.4, 0.6 ]；
% Rp                    通道纹波，如 1dB；
% Rs                    阻带纹波，如 30dB；
% 输出参数说明：
% b                     滤波器系数
% a                     滤波器系数

% Created By zxh, Date: 12-Dec-2003

Ap = (10^(Rp/20)-1)/(10^(Rp/20)+1);  % 通道纹波
As = 10^(-Rs/20);                    % 阻带衰减

% 低通滤波器
if ( Wp(1)<Ws(1) & length(Wp)==1 )
    f = [ Wp, Ws ];
    a = [ 1, 0 ];
    dev = [ Ap, As ];
end
% 带通滤波器
if ( Wp(1)>Ws(1) & length(Wp)==2 )
    f = [ Ws(1), Wp(1), Wp(2), Ws(2) ];
    a = [ 0, 1, 0 ];
    dev = [ As, Ap, As ];
end
% 高通滤波器
if ( Wp(1)>Ws(1) & length(Wp)==1 )
    f = [ Ws, Wp ];
    a = [ 0, 1 ];
    dev = [ As, Ap ];
end
% 带阻滤波器
if ( Wp(1)<Ws(1) & length(Wp)==2 )
    f = [ Wp(1), Ws(1), Ws(2), Wp(2) ];
    a = [ 1, 0, 1 ];
    dev = [ Ap, As, Ap ];
end

% 求解remez滤波器阶数
[ n, fo, ao, w ] = remezord( f, a, dev );

% remez算法
b = remez( n, fo, ao, w );  % FIR滤波器输出系数
a = 1;