function PBI = PBIGen(fc,Bw,fs,N)
% ���ɲ���Ƶ������
%     �������ã�
%     fc���������Ĺ���
%     Bandwidth�����Ŵ���
%     fs��������
%     power������
%     N�����ݳ���

%% ===================BPF �˲�����Ʋ��� ============================%

Wp = [fc-0.5*Bw  fc+0.5*Bw];       % ͨ������Ƶ��(Hz)
Ws = [fc-0.6*Bw  fc+0.6*Bw];       % �������Ƶ��(Hz)  ���ɴ�ָ��Ϊ�����10%
Rp = 3;                            % ͨ���Ʋ�(dB)
Rs = 80;                           % ���˥��(dB)

Wp = Wp/(fs/2);             % Ƶ�ʹ�һ
Ws = Ws/(fs/2);             % Ƶ�ʹ�һ

[ b, a ] = myremez( Wp, Ws, Rp, Rs );
jamFilter.Hd = dfilt.df1( b, a );                 % BPF�˲���

initialPhase = 0;                               % �ز���ʼ��λ
x = randn(N,1);    % ���ֵ��˹������
% �����ز�����
s = repmat( x, 1, length( initialPhase ) );
theta = 2*pi*fc*[1:length(s)]'/fs;
theta = repmat( theta, 1, length( initialPhase ) ) + repmat( initialPhase.', length(s), 1 );
r = s.*cos( theta );

% ���ų����˲�
jammerArray_temp1 = filter( jamFilter.Hd, r );

% ���Ź��ʹ�һ
jammerArray_temp2 = jammerArray_temp1./repmat(std(jammerArray_temp1), length(s) ,1 );

% ���ø��Ź���
PBI =  jammerArray_temp2';


return;
