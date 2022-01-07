function s = SignalGen(fc,bandwidth,Jamtype)
% �����ź�
% �������ã�  
%     Basicpara����������
%     JSR�����ű�
%     Jamtype����������
fs = 80e6;
Bw = 10e6;
Ncomp = round(2*fs/Bw);
N = Ncomp*1023;


switch (Jamtype)
    case 'S'  % ��Ƶ
        s = singletoneGen(fc,bandwidth,fs,N);
    case 'P'         % ���޸�˹
        s = PBIGen(fc,bandwidth,fs,N);
    case 'C'          % CAƥ����
        
        PRN = 5;
        numPeriod = 1;
        phase0 = 0;
        [s,am,endPhase,CAcode] = cabpsk(Ncomp,fc,fs,PRN,numPeriod,phase0);
end



end
