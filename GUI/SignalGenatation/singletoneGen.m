function singletone = singletoneGEN(fc,fs,power,N)
% ���ɵ�Ƶ������
%     �������ã�
%     fc���������Ĺ���
%     fs��������
%     power������
%     N�����ݳ���

%% ===================BPF �˲�����Ʋ��� ============================%

t=1/fs*(1:N);
ampJam  = sqrt(2*power);
singletone = ampJam*sin(2*pi*fc*t);
