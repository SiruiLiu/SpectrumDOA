function [y,am,endPhase,CAcode]=cabpsk(Ncomp,f0,Fs,PRN,numPeriod,theta)
% Ncomp��ÿ����Ԫ�ĵ���,f0���ز�,Fs�ǲ���Ƶ��,PRN�������Ǻ������������岻��������
%=========================================================================
% function������CA���BPSK�ź����ɡ�
%	����������
%	NCOMP       : number of points of each component (default: N/5)
% 	F0          : normalized frequency.              (default: 0.25)
%   PRN         : CA������
% 	Y           : signal
% 	AM          : resulting amplitude modulation     (optional).
%   numPeriod   : CA���ظ�����
% ========================================================================

%% ��ʼ������
if (nargin == 0)
    error('The number of parameters must be at least 1.');
elseif (nargin == 1)
    Ncomp=round(N/5); f0=0.25;Fs =1;
elseif (nargin == 2)
    f0=0.25;Fs =1;
end



MatlabVersion=version; MatlabVersion=str2num(MatlabVersion(1));
if (MatlabVersion==4), rand('uniform'); end


%% ����CA��

if (PRN==0)
    CAcode=ones(1,1023);
else
    CAcode = generateCAcode(PRN);
end
code = CAcode.'*ones(1,numPeriod);
code = code(:);
%% ���ɰ���
am=code*ones(1,Ncomp); %������չ��ÿ��������
am = am.';
am=am(:).';
N = length(am);
%% �����ź�
carryPhase = 2*pi*f0*(1:N)/Fs + theta;
y=am.*sin(carryPhase);
% y=sin(carryPhase);
% y=am.*exp(1i*carryPhase);
endPhase = carryPhase(end);