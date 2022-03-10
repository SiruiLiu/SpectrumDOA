function [SourceNum] = SourceEst(IQData)
    Rx = IQData*IQData'./size(IQData, 2);
    [~, D] = eig(Rx);
    Ds = diag(D);
    Factor = 0.5*(Ds(1)+Ds(2)).*eye(size(Rx));
    Rx_new = Rx - Factor;
    [V_new, D_new] = eig(Rx_new);
    G = V_new.*diag(D_new);
    W = G'.*Rx_new;
    f = mean(sqrt(abs(W)));
    delta = f - (1/size(IQData, 1))*sum(f);
    SourceNum = sum(delta > 0);
end

