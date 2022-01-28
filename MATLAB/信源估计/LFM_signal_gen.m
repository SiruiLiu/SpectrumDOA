%%--该函数产生线性调频信号，B为信号带宽，T为脉冲宽度，f0为中频，mode为线性调频模式
% mode  'R' 调频斜率为正，
%       'F' 调频斜率为负
%       'D' 调频为双边带调制
% 采样率固定为带宽的4倍
% 返回： LFM_signal： 线性调频信号
%        fs： 采样率
%        t：  对应采样时间
%        N：  采样点数
function [LFM_signal, t, N] = LFM_signal_gen(B, T, f0, mode)
    K = B/T;
    %ts = 1/fs;
    N=4096;
    
    if(mode=='D')
        t=linspace(0, T, N);
        if(0 <=t && t <T/2)
            e=exp(1i*(pi*K*t.^2 + 2*pi*f0*t));
        else
            e=exp(1i*(pi*-K*t.^2 + 2*pi*f0*t));
        end
    elseif(mode == 'R')
        t=linspace(0, T, N);
        e=exp(1i*(pi*K*t.^2 + 2*pi*f0*t));
    elseif(mode == 'F')
        t=linspace(0, T, N);
        e=exp(1i*(pi*-K*t.^2 + 2*pi*f0*t));
    end
    LFM_signal = e;
    
end

