function [y,am,endPhase,CAcode]=cabpsk(Ncomp,f0,Fs,PRN,numPeriod,theta)
% Ncomp是每个码元的点数,f0是载波,Fs是采样频率,PRN是温卫星号其他参数含义不懂在问我
%=========================================================================
% function：具有CA码的BPSK信号生成。
%	参数描述：
%	NCOMP       : number of points of each component (default: N/5)
% 	F0          : normalized frequency.              (default: 0.25)
%   PRN         : CA码的序号
% 	Y           : signal
% 	AM          : resulting amplitude modulation     (optional).
%   numPeriod   : CA码重复周期
% ========================================================================

%% 初始化设置
if (nargin == 0)
    error('The number of parameters must be at least 1.');
elseif (nargin == 1)
    Ncomp=round(N/5); f0=0.25;Fs =1;
elseif (nargin == 2)
    f0=0.25;Fs =1;
end



MatlabVersion=version; MatlabVersion=str2num(MatlabVersion(1));
if (MatlabVersion==4), rand('uniform'); end


%% 生成CA码

if (PRN==0)
    CAcode=ones(1,1023);
else
    CAcode = generateCAcode(PRN);
end
code = CAcode.'*ones(1,numPeriod);
code = code(:);
%% 生成包络
am=code*ones(1,Ncomp); %将码扩展到每个采样点
am = am.';
am=am(:).';
N = length(am);
%% 生成信号
carryPhase = 2*pi*f0*(1:N)/Fs + theta;
y=am.*sin(carryPhase);
% y=sin(carryPhase);
% y=am.*exp(1i*carryPhase);
endPhase = carryPhase(end);