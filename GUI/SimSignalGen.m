function [SimSignal] = SimSignalGen(SigNum, DOA,SigInfo,Snapshot, ChannelSNR, ArrayPos, Lambda)
    DOA(all(DOA == 0, 2),:) = [];
    SigInfo(all(SigInfo==0, 2), :) = [];
    Lambda(Lambda == 0) = [];
    SNRdB = SigInfo(:, 3);
    DOA = DOA*pi/180;

%     Sig = zeros(SigNum, 16368);
%     for i = 1:SigNum
%         Sig(i,:) = SignalGen(SigInfo(i, 1), SigInfo(i, 2), 'C');
%     end
    
    Sig = zeros(SigNum, 16368);
    for i = 1:SigNum
        Sig(i,:) = SignalGen(SigInfo(i, 1), SigInfo(i, 2), 'P');
        Sig(i,:) = awgn(Sig(i,:), SNRdB(i), 'measured');
    end
    Sig = Sig(:,1:Snapshot);
    [xs, xn] = arraydata(Sig, ChannelSNR, DOA, ArrayPos, Lambda);
    SimSignal = xs+xn;
end

