function s = SignalGen(fc,bandwidth,Jamtype)
% 生成信号
% 参数设置：  
%     Basicpara：基本参数
%     JSR：干信比
%     Jamtype：干扰类型
fs = 80e6;
Bw = 10e6;
Ncomp = round(2*fs/Bw);
N = Ncomp*1023;


switch (Jamtype)
    case 'S'  % 单频
        s = singletoneGen(fc,bandwidth,fs,N);
    case 'P'         % 带限高斯
        s = PBIGen(fc,bandwidth,fs,N);
    case 'C'          % CA匹配谱
        
        PRN = 5;
        numPeriod = 1;
        phase0 = 0;
        [s,am,endPhase,CAcode] = cabpsk(Ncomp,fc,fs,PRN,numPeriod,phase0);
end



end
