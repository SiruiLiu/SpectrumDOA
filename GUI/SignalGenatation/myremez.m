function [ b, a ] = myremez( Wp, Ws, Rp, Rs )
% ����parks-mcclellan��������FIR�˲����Ľ���������Remez�㷨��������FIR�˲������
% �������˵����
% Wp                    ͨ������Ƶ�ʣ���[ 0.3, 0.7 ]�����巽����matlab��ͬ
% Ws                    �������Ƶ�ʣ���[ 0.4, 0.6 ]��
% Rp                    ͨ���Ʋ����� 1dB��
% Rs                    ����Ʋ����� 30dB��
% �������˵����
% b                     �˲���ϵ��
% a                     �˲���ϵ��

% Created By zxh, Date: 12-Dec-2003

Ap = (10^(Rp/20)-1)/(10^(Rp/20)+1);  % ͨ���Ʋ�
As = 10^(-Rs/20);                    % ���˥��

% ��ͨ�˲���
if ( Wp(1)<Ws(1) & length(Wp)==1 )
    f = [ Wp, Ws ];
    a = [ 1, 0 ];
    dev = [ Ap, As ];
end
% ��ͨ�˲���
if ( Wp(1)>Ws(1) & length(Wp)==2 )
    f = [ Ws(1), Wp(1), Wp(2), Ws(2) ];
    a = [ 0, 1, 0 ];
    dev = [ As, Ap, As ];
end
% ��ͨ�˲���
if ( Wp(1)>Ws(1) & length(Wp)==1 )
    f = [ Ws, Wp ];
    a = [ 0, 1 ];
    dev = [ As, Ap ];
end
% �����˲���
if ( Wp(1)<Ws(1) & length(Wp)==2 )
    f = [ Wp(1), Ws(1), Ws(2), Wp(2) ];
    a = [ 1, 0, 1 ];
    dev = [ Ap, As, Ap ];
end

% ���remez�˲�������
[ n, fo, ao, w ] = remezord( f, a, dev );

% remez�㷨
b = remez( n, fo, ao, w );  % FIR�˲������ϵ��
a = 1;