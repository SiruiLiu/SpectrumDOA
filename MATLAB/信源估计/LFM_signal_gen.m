%%--�ú����������Ե�Ƶ�źţ�BΪ�źŴ���TΪ�����ȣ�f0Ϊ��Ƶ��modeΪ���Ե�Ƶģʽ
% mode  'R' ��Ƶб��Ϊ����
%       'F' ��Ƶб��Ϊ��
%       'D' ��ƵΪ˫�ߴ�����
% �����ʹ̶�Ϊ�����4��
% ���أ� LFM_signal�� ���Ե�Ƶ�ź�
%        fs�� ������
%        t��  ��Ӧ����ʱ��
%        N��  ��������
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

